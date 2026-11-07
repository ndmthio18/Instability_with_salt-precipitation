"""Package source-only SUNDIALS from its pinned upstream release; no builds here."""
import hashlib
import pathlib
import tarfile
import sys

root = pathlib.Path(sys.argv[1])
archive = root / 'sundials-6.7.0.tar.gz'
expected = 'd434d00a68d0f18dc2c8f4c306e4becf62ac6d37a259253ec9ea616152babba4'
assert hashlib.sha256(archive.read_bytes()).hexdigest() == expected
allowed = {'CMakeLists.txt', 'cmake', 'include', 'src', 'LICENSE', 'NOTICE', 'README.md'}
with tarfile.open(archive) as source, tarfile.open(root / 'sundials-source.tar.gz', 'w:gz') as dest:
    for member in source.getmembers():
        bits = pathlib.PurePosixPath(member.name).parts
        if len(bits) == 1 or bits[1] in allowed:
            dest.addfile(member, source.extractfile(member) if member.isfile() else None)
print(hashlib.sha256((root / 'sundials-source.tar.gz').read_bytes()).hexdigest())
