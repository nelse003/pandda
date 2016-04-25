import os, sys, shutil

#:::::::::::::::::::::::::::::::::#
# ############################### #
# ### Path Utility Functions  ### #
# ############################### #
#:::::::::::::::::::::::::::::::::#

def easy_directory(directory):
    """Checks a directory exists and creates it if not"""
    if not os.path.exists(directory):
        os.mkdir(directory)
    return directory

def rel_symlink(orig, link):
    """Make a relative symlink from link to orig"""
    assert os.path.exists(orig), 'FILE DOES NOT EXIST: {!s}'.format(orig)
    assert not os.path.exists(link), 'LINK ALREADY EXISTS: {!s}'.format(link)
    orig = os.path.abspath(orig)
    link = os.path.abspath(link)
    assert not link.endswith('/'), 'LINK CANNOT END WITH /'
    os.symlink(os.path.relpath(orig, start=os.path.dirname(link)), link)

def delete_temporary_directory(tempdir):
    """Delete A Temporary Directory"""

    # Remove Temporary Files
    if tempdir.startswith('/tmp/'):
        shutil.rmtree(tempdir)
        return 0
    else:
        FlagError("Will Not Delete This Directory - Not in '/tmp' - Remove Manually : "+tempdir)
        return 1

def list_directory_contents(directory, templatename=None, templatestyle=None):
    """List Folders in a Directory"""

    # List the contents of the directory
    contents = os.listdir(directory)
    # Decide which are directories and which to keep
    if   templatename and templatestyle=='End':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and dir.endswith(templatename)]
    elif templatename and templatestyle=='Start':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and dir.startswith(templatename)]
    elif templatename and templatestyle=='Contains':
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir)) and templatename in dir]
    else:
        subdirs = [dir for dir in contents if os.isdir(os.path.join(directory,dir))]
    # Construct Full Paths
    subdirs = [os.path.join(directory,dir) for dir in subdirs]

    return subdirs
