# Contributing

Please open a GitHub issue for a reproducible bug or a proposed behavioral change.
Include the RaMA-G version, command, platform, compiler and relevant log. A small
artificial input is preferable to a private or very large genome. Do not post
credentials, private paths or sequence data without permission.

Build with the `testing` preset and run CTest with `--no-tests=error`. Add a small
independent correctness check for coordinate, CIGAR, score or ownership changes.
Keep production behavior changes separate from file moves and formatting.
Document changes to supported command-line or installed library interfaces.

Retain third-party license notices and source provenance. Do not edit vendored
upstream snapshots merely to change this project's layout. Public CI uses only
synthetic inputs; large research datasets and experiment archives are not part
of a contribution. See the [developer guide](docs/developer-guide.md).
