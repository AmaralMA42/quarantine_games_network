"""Attach the retained scientific data bundle without replacing repository sources."""
from pathlib import Path, PurePosixPath
import csv
import hashlib
import io
import json
import stat
import zipfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def contained_path(root, name):
    """Keep every destination below the selected repository root."""
    path = PurePosixPath(name)
    if path.is_absolute() or '\\' in name or ':' in name or any(x in {'..', '.'} for x in path.parts):
        raise ValueError('Invalid bundle member path')
    target = root.joinpath(*path.parts)
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError('Bundle path leaves the repository')
    return target


def _prepare_base(archive, check_only=False):
    root = ROOT.resolve()
    archive = Path(archive).resolve()
    specification = json.loads((root/'DATA_BUNDLE.json').read_text(encoding='utf-8'))
    if digest(archive) != specification['sha256']:
        raise ValueError('The archive checksum differs from DATA_BUNDLE.json')
    tracked = {row['path'] for row in csv.DictReader((root/'SOURCE_MANIFEST_SHA256.csv').open(encoding='utf-8'))}
    tracked.update({'SOURCE_MANIFEST_SHA256.csv', '.gitignore'})
    prefix = specification['zip_root'] + '/'
    # Verification precedes any filesystem change. Existing scientific inputs
    # are preserved; mismatches stop the operation rather than overwrite them.
    with zipfile.ZipFile(archive) as bundle:
        infos = bundle.infolist()
        if len({x.filename for x in infos}) != len(infos):
            raise ValueError('Duplicate archive member')
        for info in infos:
            if not info.filename.startswith(prefix):
                raise ValueError('Unexpected archive root')
            if stat.S_ISLNK(info.external_attr >> 16):
                raise ValueError('Symlink member is not supported')
            contained_path(root, info.filename[len(prefix):])
        rows = list(csv.DictReader(io.StringIO(bundle.read(prefix+'MANIFEST_SHA256.csv').decode('utf-8-sig'))))
        inputs = list(csv.DictReader(io.StringIO(bundle.read(prefix+'provenance/SELECTED_INPUTS.csv').decode('utf-8-sig'))))
        if len(rows) != specification['manifest_files'] or len(inputs) != specification['original_bindings']:
            raise ValueError('Unexpected scientific manifest size')
        paths = [row['path'] for row in rows]
        if len(set(paths)) != len(paths):
            raise ValueError('Duplicate manifest path')
        if set(paths) | {'MANIFEST_SHA256.csv'} != {info.filename[len(prefix):] for info in infos}:
            raise ValueError('Manifest and archive members differ')
        generated_prefixes = ('reproduction/main/manuscript/figures/', 'reproduction/si/supplement/figures/')
        immutable = {row['path'] for row in inputs if not row['path'].startswith(generated_prefixes)}
        planned = []
        protected, existing, kept_generated = [], [], []
        for row in rows:
            name = row['path']
            target = contained_path(root, name)
            data = bundle.read(prefix + name)
            if len(data) != int(row['bytes']) or hashlib.sha256(data).hexdigest() != row['sha256']:
                raise ValueError('Scientific bundle manifest mismatch: ' + name)
            if name in tracked:
                protected.append(name)
                continue
            if target.exists():
                if not target.is_file():
                    raise FileExistsError('A directory occupies an expected file: ' + name)
                if digest(target) == row['sha256']:
                    existing.append(name)
                elif name in immutable:
                    raise FileExistsError('An existing scientific input differs: ' + name)
                else:
                    kept_generated.append(name)
            else:
                planned.append(row)
        result = dict(status='CHECKED' if check_only else 'PREPARED',
                      archive_sha256=specification['sha256'],
                      manifest_files_checked=len(rows), original_bindings=len(inputs),
                      protected_repository_files=len(protected),
                      unchanged_existing_files=len(existing),
                      preserved_generated_files=len(kept_generated),
                      files_to_add=len(planned), files_added=0,
                      source_files_overwritten=0, new_epidemic_runs=0,
                      network_download=False)
        if check_only:
            return result
        for row in planned:
            target = contained_path(root, row['path'])
            target.parent.mkdir(parents=True, exist_ok=True)
            with target.open('xb') as stream:
                stream.write(bundle.read(prefix + row['path']))
            if digest(target) != row['sha256']:
                raise ValueError('Newly copied file failed verification: ' + row['path'])
        result['files_added'] = len(planned)
    (root/'qa').mkdir(exist_ok=True)
    (root/'qa/DATA_ATTACHMENT.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    return result


def _prepare_update(archive, check_only=False):
    root=ROOT.resolve();archive=Path(archive).resolve()
    specification=json.loads((root/'ARTICLE_UPDATE.json').read_text('utf-8'))
    if digest(archive)!=specification['sha256']:raise ValueError('Current update archive checksum mismatch')
    prefix=specification['zip_root']+'/'
    manifest_name=specification['manifest_path']
    with zipfile.ZipFile(archive) as bundle:
        infos=bundle.infolist()
        if len({x.filename for x in infos})!=len(infos):raise ValueError('Duplicate update member')
        for info in infos:
            if not info.filename.startswith(prefix):raise ValueError('Unexpected update root')
            name=info.filename[len(prefix):]
            if not name.startswith('article_update/'):raise ValueError('Update may only write inside article_update')
            if stat.S_ISLNK(info.external_attr>>16):raise ValueError('Update symlink is not supported')
            contained_path(root,name)
        manifest=bundle.read(prefix+manifest_name)
        if hashlib.sha256(manifest).hexdigest()!=specification['manifest_sha256']:raise ValueError('Update manifest checksum mismatch')
        rows=list(csv.DictReader(io.StringIO(manifest.decode('utf-8-sig'))))
        paths=[x['path'] for x in rows]
        if len(paths)!=len(set(paths)) or len(rows)!=specification['manifest_files']:raise ValueError('Unexpected update manifest')
        if set(paths)|{manifest_name}!={x.filename[len(prefix):] for x in infos}:raise ValueError('Update manifest/member mismatch')
        planned=[]
        for row in rows+[{'path':manifest_name,'bytes':len(manifest),'sha256':specification['manifest_sha256']}]:
            name=row['path'];target=contained_path(root,name);data=bundle.read(prefix+name)
            if len(data)!=int(row['bytes']) or hashlib.sha256(data).hexdigest()!=row['sha256']:raise ValueError('Update member mismatch: '+name)
            if target.exists():
                if not target.is_file() or digest(target)!=row['sha256']:raise FileExistsError('Current input differs: '+name)
            else:planned.append(row)
        if not check_only:
            for row in planned:
                target=contained_path(root,row['path']);target.parent.mkdir(parents=True,exist_ok=True)
                with target.open('xb') as stream:stream.write(bundle.read(prefix+row['path']))
        return dict(status='CHECKED' if check_only else 'PREPARED',archive_sha256=specification['sha256'],
            manifest_files_checked=len(rows),files_added=0 if check_only else len(planned),source_files_overwritten=0)


def prepare(archive, check_only=False, update_archive=None):
    # Verify both inputs before adding files. Historical bindings stay intact;
    # the current figure selection has its own namespace and manifest.
    if update_archive is not None:
        _prepare_base(archive,True)
        _prepare_update(update_archive,True)
    result=_prepare_base(archive,check_only)
    if update_archive is not None:result['current_figure_update']=_prepare_update(update_archive,check_only)
    else:result['current_figure_update']='NOT_ATTACHED; provide --update-archive for the current article'
    return result
