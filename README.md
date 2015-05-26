## Introduction

Hi there!  You've just downloaded yt, an analysis tool for scientific
datasets, generated on a variety of data platforms.  It's written in 
python and heavily leverages NumPy, Matplotlib, SymPy and Cython for a variety
of tasks.

Full documentation and a user community can be found at:

http://yt-project.org/

http://yt-project.org/doc/

If you have used Python before, and are comfortable with installing packages,
you should find the setup.py script fairly straightforward: simply execute
"python setup.py install".

If you would rather a more automated installation, you can use the script
doc/install_script.sh .  You will have to set the destination directory, and
there are options available, but it should be straightforward.

For more information on installation, what to do if you run into problems, or 
ways to help development, please visit our website.

Enjoy!

## Required Packages

If you're using a Linux distro, you will need to have these packages installed. Use your distro's package distro to install:

 - hdf5
 - zeromq
 - sqlite
 - mercurial (this is what we use for version control)

On [Arch Linux](https://archlinux.org) these can be installed with `pacman`:

```
$ pacman -S hdf5 zeromq sqlite mercurial
```

On [Gentoo](https://gentoo.org) with `emerge`:

```
$ emerge -a hdf5 zeromq sqlite mercurial
```

On [Debian](http://debian.org)/[Ubuntu](http://ubuntu.com):

```
$ apt-get install build-essential hd5 zeromq sqlite mercurial
```

## Setup and Installation with pyenv

[`pyenv`](https://github.com/yyuu/pyenv) is a python version manager. See [here](https://github.com/yyuu/pyenv/blob/master/README.md) for info on installation and setup.

I would also recommend installing [`pyenv-pip-rehash`](https://github.com/yyuu/pyenv-pip-rehash) to prevent having to type `pyenv rehash` every time you install something with `pip` or a python version.

Then clone this repo locally with mercurial (`hg`):

```
$ hg clone https://bitbucket.org/astrohckr/yt
```

Change into the `yt` directory:

```
$ cd /path/to/yt
```

ex. I cloned this repo to my `~/code` folder:

```
$ cd ~/code/yt
```

When you change into the `yt` directory, `pyenv` will scan it for a `.python-version` file that specifies what version of python the project uses. For our use case, this is `2.7.9`. 

To install that version of python, you simply have to run:

```
$ pyenv install
```

Once that has completed, you can now use `pip` to install the required python packages with the `requirements.txt` file in the root `yt` path:

```
$ pip install -r requirements.txt
```

Then finally, to setup `yt`:

```
$ python setup.py install
```



