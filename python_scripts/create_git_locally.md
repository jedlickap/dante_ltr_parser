## Create git repository 
### - locally and on GitHub and push all scripts there

```bash

# init
(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git init
(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ vi .gitignore 
$ cat .gitignore 
# Ignore Python compiled files and caches
__pycache__/
*.pyc

# Ignore virtual environment folder
dante_ltr_parser/

# Ignore OS-specific files
.DS_Store

# Create project on github
# Go to GitHub and log into your account.
# Click the + icon in the top right corner, and select New repository.
# Give your repository a name, such as dante_ltr_parser, and add an optional description.
# Choose whether you want the repository to be public or private.
# Do not initialize the repository with a README (since you already have files).
# Click Create repository.

(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git add .
(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ echo "# dante_ltr_parser" >> README.md

(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git add README.md 

(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git commit -m "Initial commit of dante_ltr_parser project"
[master 08f8083] Initial commit of dante_ltr_parser project
 1 file changed, 1 insertion(+)
 create mode 100644 README.md
(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git branch -M main
(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git remote add origin https://github.com/jedlickap/dante_ltr_parser.git
(dante_ltr_parser) pavel@debian11:~/Documents/dante_ltr_parser$ git push -u origin main
Enumerating objects: 25, done.
Counting objects: 100% (25/25), done.
Delta compression using up to 8 threads
Compressing objects: 100% (20/20), done.
Writing objects: 100% (25/25), 12.23 KiB | 6.12 MiB/s, done.
Total 25 (delta 4), reused 0 (delta 0), pack-reused 0
remote: Resolving deltas: 100% (4/4), done.
To https://github.com/jedlickap/dante_ltr_parser.git
 * [new branch]      main -> main
Branch 'main' set up to track remote branch 'main' from 'origin'.
```

### Update git repository

```
# Check the status of the repository:
$ git status

# Stage all changes (or specific files):

$ git add .

# or for specific files:
$ git add path/to/modified_file.py

# Commit the changes:
$ git commit -m "Your descriptive commit message"

# Push the changes to GitHub:
$ git push origin main
```
