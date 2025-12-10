"""
Set your local repository path for saving/reading data.
Replace the placeholder below with an absolute path to this clone, e.g.
projectdir = "C:/Users/yourname/projects/quantum-mcmc/"
"""
from pathlib import Path

# TODO: update this to your local absolute path if you want to hard-code it
projectdir = "C:/path/to/quantum-mcmc/"

if projectdir.startswith("C:/path/to"):
    repo_root = Path(__file__).resolve().parent.parent
    print(
        "Please set src/path.py:projectdir to your local absolute path "
        f"(default suggestion: {repo_root}/)"
    )
