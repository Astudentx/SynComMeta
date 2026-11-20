"""Plan paired-read analyses and export self-contained Bash jobs."""
import argparse
import csv
import fcntl
import hashlib
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from .modules import alignment, coverage, m01, m02, m03, m04, m05, m06

HERE = Path(__file__).resolve().parent
MODULES = {
    'M01': (m01, 'fastp 配对 reads 质控'),
    'M02': (m02, 'mOTUs 分类及丰度'),
    'M03': (m03, '参考菌株及近缘组覆盖'),
    'M04': (m04, 'CDS 片段计数与已有功能注释汇总'),
    'M05': (m05, '群体 SNV 与可调用区间'),
    'M06': (m06, '基因 PAV 与相对拷贝数'),
}
TOOLS = {'M01': ['fastp'], 'M02': ['motus', 'bwa'], 'M03': ['bwa', 'samtools'],
         'M04': ['bwa', 'samtools', 'featureCounts'], 'M05': ['bcftools'], 'M06': []}


def table(path):
    with path.open(newline='') as handle:
        yield from csv.DictReader(handle, delimiter='\t')


def pairs(path):
    result = {}
    for line in path.read_text().splitlines():
        if line and not line.startswith('#'):
            key, value = line.split('\t')
            if not value or key in result:
                raise ValueError(f'Invalid entry in {path}: {key}')
            result[key] = value
    return result


def require(path):
    if not path.is_file() or not path.stat().st_size or not os.access(path, os.R_OK):
        raise ValueError(f'Missing or empty input: {path}')
    return path


def samples(path, clean):
    reads = defaultdict(list)
    seen = set()
    for row in table(path):
        sid = row['sample_id']
        lane = row.get('lane_id') or '.'
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or not re.fullmatch(r'[\w.-]+', lane):
            raise ValueError(f'Invalid sample/lane: {sid}/{lane}')
        if (sid, lane) in seen or (reads[sid] and (clean or lane == '.' or (sid,'.') in seen)):
            raise ValueError(f'Duplicate sample or ambiguous lanes: {sid}')
        mates = [require((path.parent / row[key]).resolve()) for key in ('r1','r2')]
        if mates[0].samefile(mates[1]):
            raise ValueError(f'Identical mates: {sid}')
        seen.add((sid,lane))
        reads[sid].append((lane,*mates))
    if not reads:
        raise ValueError('Empty sample table')
    return reads


def database(path, selected):
    ref = pairs(path / 'manifest.tsv')
    version = ref.pop('db_version')
    if not re.fullmatch(r'[\w.-]+', version):
        raise ValueError('Invalid database version')
    ref = {k: str((path / v).resolve()) for k,v in ref.items()}
    if 'M02' in selected:
        require(Path(ref['motus_db'])/'db_mOTU_versions')
    if selected & {'M03','M04'}:
        for ext in ('', '.amb','.ann','.bwt','.pac','.sa','.fai'):
            require(Path(ref['reference']+ext))
        rows = list(table(require(Path(ref['contig_map']))))
        lengths = {r['contig_id']: int(r['length_bp']) for r in rows}
        with Path(ref['reference']+'.fai').open() as handle:
            index = [line.split('\t')[:2] for line in handle]
        if len(rows)!=len(lengths) or len(index)!=len(rows) or lengths!={k:int(v) for k,v in index}:
            raise ValueError('Reference and contig map disagree')
        if 'strains' in ref:
            records = list(table(require(Path(ref['strains']))))
            ids = {r['ID'] for r in records if r['strain_id']}
            if len(ids)!=len(records) or ids!={r['genome_id'] for r in rows}:
                raise ValueError('Strain metadata and contig map disagree')
    for mod, keys in [('M04',['genes_saf','annotations']),('M06',['genes_bed'])]:
        if mod in selected:
            for key in keys:
                require(Path(ref[key]))
    return dict(ref, db_version=version)


def generate(plan, args, reads, selected, tools, ref):
    quote = lambda value: shlex.quote(str(value))
    out = args.outdir.resolve()
    runtime_keys = ("db_version", "reference", "contig_map", "motus_db", "genes_saf", "genes_bed", "annotations")
    settings = dict({k: ref[k] for k in runtime_keys if k in ref}, **tools, out=out, threads=args.threads, clean=int(args.clean_reads),
                    min_genome_depth=args.min_genome_depth, min_genome_breadth=args.min_genome_breadth)
    settings.update({k:int(k in selected) for k in MODULES})
    (plan/'settings.sh').write_text(''.join(f'{k}={quote(v)}\n' for k,v in settings.items()
                                           if re.fullmatch(r'[a-zA-Z_][a-zA-Z0-9_]*',k)))
    shutil.copyfile(args.samples, plan/'samples.input.tsv')
    shutil.copyfile(args.db/'manifest.tsv', plan/'database.tsv')
    if ref.get('strains'):
        shutil.copyfile(ref['strains'], plan/'strains.tsv')
    order = [MODULES[k][0] for k in ('M01','M02') if k in selected]
    if selected & {'M03','M04'}:
        order.append(alignment)
    order.extend(MODULES[k][0] for k in ('M03','M04') if k in selected)
    if selected & {'M05','M06'}:
        order.append(coverage)
    order.extend(MODULES[k][0] for k in ('M05','M06') if k in selected)
    commands = [module.COMMANDS for module in order]
    if 'M03' in selected and coverage not in order:
        commands.append(coverage.COMMANDS)
    commands.append((HERE/'tables.bash').read_text())
    (plan/'commands.sh').write_text('\n'.join(commands))
    runtime = (HERE/'runtime.bash').read_text()
    summary = (HERE/'summary.bash').read_text()
    digest = hashlib.sha256()
    for file in sorted(p for p in HERE.rglob('*') if p.suffix in {'.py','.bash'}):
        digest.update(file.read_bytes())
    (plan/'method.sha256').write_text(digest.hexdigest()+'\n')
    with (plan/'reads.tsv').open('w') as manifest:
        for sid in sorted(reads):
            lanes = reads[sid]
            for row in lanes:
                print(sid,*row,sep='\t',file=manifest)
            job = ['sample='+quote(sid), 'qc="$out/qc/$sample"']
            for i, key in enumerate(('r1','r2'),1):
                job.append(key+'=('+' '.join(quote(row[i]) for row in lanes)+')')
            for module in order:
                if module in (m01,m02,alignment):
                    job.append('reads_signature=$(identity "${r1[@]}" "${r2[@]}")')
                job.append(module.STEP)
            name = f'sample_{sid}.sh'
            (plan/name).write_text('\n'.join(job)+'\n')
            runtime += f'\nsource "$plan/{name}"\n'
    (plan/'run_all.sh').write_text(runtime+'\n'+summary)
    (plan/'run_all.sh').chmod(0o755)


def execute(args):
    if args.action == 'modules':
        for key, (_, label) in MODULES.items():
            print(key, label, sep='\t')
        return
    configured = pairs(args.tools)
    if set(configured)-{t for names in TOOLS.values() for t in names}:
        raise ValueError('Unknown tool in configuration')
    selected = set(MODULES) if args.action=='doctor' else set(args.modules.split(','))
    if selected-set(MODULES):
        raise ValueError('Unknown module')
    if args.action!='doctor':
        if args.clean_reads and 'M01' in selected:
            raise ValueError('M01 conflicts with --clean-reads')
        if not args.clean_reads:
            selected.add('M01')
        if selected & {'M05','M06'}:
            selected.add('M03')
    tools = {}
    for key in sorted({t for m in selected for t in TOOLS[m]}):
        value = configured.get(key,key)
        if '/' in value:
            value = str((args.tools.parent/value).absolute())
        path = shutil.which(value)
        if not path:
            raise ValueError(f'Missing executable: {key} ({value})')
        tools[key] = path
        if args.action=='doctor':
            print(key,path,sep='\t')
    if args.action=='doctor':
        return
    if args.threads<1 or not (0<args.min_genome_depth<float('inf') and 0<args.min_genome_breadth<=1):
        raise ValueError('Invalid thread count or coverage threshold')
    reads = samples(args.samples, args.clean_reads)
    ref = database(args.db, selected)
    out = args.outdir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    with (out/'.syncom.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        plans = out/'.plans'
        plans.mkdir(exist_ok=True)
        plan = Path(tempfile.mkdtemp(prefix='python-', dir=plans))
        generate(plan,args,reads,selected,tools,ref)
        for name in ('scripts','shell'):
            link = plan.with_suffix('.'+name)
            link.symlink_to(plan, target_is_directory=True)
            link.replace(out/name)
        print(out/'shell/run_all.sh',flush=True)
        if args.action=='run':
            # Bash uses fd 9 for its output lock, shared with the planner.
            os.dup2(lock.fileno(),9,inheritable=True)
            os.set_inheritable(9,True)
            os.execvp('bash',['bash',str(plan/'run_all.sh'),'--lock-held'])


def main():
    if sys.argv[1:2]==['merge']:
        sys.exit(subprocess.call(['bash',str(HERE/'merge.bash'),*sys.argv[2:]]))
    parser=argparse.ArgumentParser(description='SynCom-MetaFlow Python workflow (merge --outdir DIR BATCH BATCH ...)')
    sub=parser.add_subparsers(dest='action',required=True)
    sub.add_parser('modules')
    for name in ('run','plan','doctor'):
        cli=sub.add_parser(name)
        cli.add_argument('--tools',type=Path,default=HERE.parent/'config/tools.tsv')
        if name!='doctor':
            for option in ('samples','db','outdir'):
                cli.add_argument('--'+option,type=Path,required=True)
            cli.add_argument('--modules',required=True)
            cli.add_argument('--threads',type=int,default=2)
            cli.add_argument('--min-genome-depth',type=float,default=5)
            cli.add_argument('--min-genome-breadth',type=float,default=0.5)
            cli.add_argument('--clean-reads',action='store_true')
    try:
        execute(parser.parse_args())
    except (ValueError, KeyError, OSError) as exc:
        parser.exit(1,f'Error: {exc}\n')
