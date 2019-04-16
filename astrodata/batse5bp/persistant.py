"""
Utility classes for simple object persistance.

More sophisticated alternatives include atpy, h5py, and PyTables.  This
simple approach is adopted for now to reduce dependencies.

*** MODULE NOT YET USED ***

Based on ioutils.py from the Inference package.

Created 2012-10-22 by Tom Loredo
2019:  Converted to Python 3
"""

# TODO:  Import backup

import pickle as pkl
from numpy import ndarray, savez, load

__all__ = ['ArrayStore', 'AttrStore']


class ArrayStore:
    """
    A container class that provides lazy persistance of array instance
    attributes via NumPy's .npz archives.

    This is doubly lazy:  arrays are not read into the namespace until
    actually requested, and they are not saved to the archive until the
    ArrayStore is explicitly saved.
    """

    # *** Is loading really lazy, i.e., does NumPy's load() read in array data
    # before actual access via dict lookup?

    # TODO:  Perhaps have a __del__ method that saves on deletion?
    # It would have to handle fname=None gracefully.
    # Perhaps an init arg should specify whether to save on deletion.

    # TODO:  Support creating with a name to a non-existing file;
    # subsequent "save()" will use the stored name.

    def __init__(self, fname=None):
        """
        Prepare to load arrays from storage if a file name is provided;
        otherwise support saving of arrays assigned as attributes.
        """
        # All internal attributes start with '_' to avoid __setattr__
        # array filtering.
        if fname:
            if not fname.endswith('.npz'):
                fname = fname + '.npz'
                self._npz = load(fname)
        else:
            self._npz = None
        self._fname = fname
        self._archived = {}  # arrays pulled from the archive
        self._new = {}  # arrays to be archived

    def __getattr__(self, name):
        """
        Catch references to array attributes that have not yet been
        loaded from the archive.
        """
        # This is called only if name is not already in the instance dict.
        if self._npz is None:  # no archive to grab attribute from
            raise AttributeError(name)
        else:  # get value from archive and keep a reference to it
            try:
                value = self._npz[name]
                # Set the attribute directly so __setattr__ won't add it
                # to self.new.
                object.__setattr__(self, name, value)
                self._archived[name] = value
                return value
            except:
                raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        Catch assignments to new attributes, marking them for saving when
        the store is next saved.

        Names starting with '_' have their corresponding attributes set
        without marking; such names are intended for internal use only,
        to keep track of the state of the store.
        """
        # *** Should this prevent over-writing existing names, either
        # directly assigned or yet to be loaded from the archive?
        # Right now we allow reassignments; saving will save the
        # reassigned value.
        if name.startswith('_'):
            object.__setattr__(self, name, value)
        elif isinstance(value, ndarray):
            self._new[name] = value
            object.__setattr__(self, name, value)
        else:
            raise ValueError('Only ndarray objects may be stored!')

    def save(self, fname=None):
        """
        Save array attributes to a NumPy .npz archive.  If a name is
        provided, it is used (adding '.npz' if needed); otherwise it is
        presumed this store was created from an existing archive, and
        that archive's name is used.

        In either case, any existing version is backed up before the
        new one is created, with up to two levels of backups (with '%'
        and '%%' suffixes).
        """
        # *** Perhaps this should be called by __del__; or this possibility
        # could be set by a flag on init???

        # If no name given, use the name provided at creation.
        if fname is None:
            if self._fname is None:
                raise ValueError('Need a file name for saving!')
            fname = self._fname
        else:
            if not fname.endswith('.npz'):
                fname = fname + '.npz'
        # Gather arrays to store; note new versions supersede old ones.
        arrays = {}
        for name, val in self._archived.items():
            arrays[name] = val
        # Don't forget to get any archived values not already accessed.
        if self._npz:
            for name in self._npz.files:
                if name not in self._archived:
                    arrays[name] = self._npz[name]
        for name, val in self._new.items():
            arrays[name] = val
        # Backup any existing file, and save the arrays to a new file.
        backup(fname)
        savez(fname, **arrays)

    def contents(self):
        """
        List the names of arrays accessible by this store, including
        both previously archived arrays and newly defined arrays.
        """
        names = []
        if self._npz:
            names.extend(self._npz.files)  # archive file names will be array names
        names.extend(list(self._new.values()))
        return names


class AttrStore:
    """
    A container class that provides lazy persistance of instance attributes
    as a Python pickle.

    This is doubly lazy:  attributes are not placed in the namespace until
    actually requested (but they are all read from the file), and they are not
    saved to the archive until the AttrStore is explicitly saved.

    If only NumPy arrays are to be stored, ArrayStore may be more efficient.
    """

    # *** Is there any virtue to laziness in loading the namespace, since
    # we aren't unpacking a NumPy npz file archive?

    # TODO:  Perhaps have a __del__ method that saves on deletion?
    # It would have to handle fname=None gracefully.
    # Perhaps an init arg should specify whether to save on deletion.

    # TODO:  Support creating with a name to a non-existing file;
    # subsequent "save()" will use the stored name.

    # TODO:  Add kwd args to init for setting attributes to store.

    def __init__(self, fname=None, **kwds):
        """
        Prepare to load attributes from storage if a file name is provided;
        otherwise support saving of arrays assigned as attributes.
        """
        # All internal attributes start with '_' to avoid __setattr__
        # array filtering.
        if fname:
            if not fname.endswith('.pkl'):
                fname = fname + '.pkl'
                ifile = open(fname, 'rb')
                self._archive = pkl.load(ifile)
                ifile.close()
        else:
            self._archive = None
        self._fname = fname
        self._pulled = {}  # attributes pulled from the archive
        self._new = {}  # attributes to be archived
        if kwds:
            for key, value in kwds.items():
                setattr(self, key, value)

    def __getattr__(self, name):
        """
        Catch references to attributes that have not yet been loaded from
        the store.
        """
        # This is called only if name is not already in the instance dict.
        if self._archive is None:  # no archive to grab attribute from
            raise AttributeError(name)
        else:  # get value from archive and keep a reference to it
            try:
                value = self._archive[name]
                # Set the attribute directly so __setattr__ won't add it
                # to self.new.
                object.__setattr__(self, name, value)
                self._pulled[name] = value
                return value
            except:
                raise AttributeError(name)

    def __setattr__(self, name, value):
        """
        Catch assignments to new attributes, marking them for saving when
        the store is next saved.

        Names starting with '_' have their corresponding attributes set
        without marking; such names are intended for internal use only,
        to keep track of the state of the store.
        """
        # *** Should this prevent over-writing existing names, either
        # directly assigned or yet to be loaded from the archive?
        # Right now we allow reassignments; saving will save the
        # reassigned value.
        if name.startswith('_'):
            object.__setattr__(self, name, value)
        else:
            self._new[name] = value
            object.__setattr__(self, name, value)

    def save(self, fname=None):
        """
        Save attributes to a high-protocal pickle.  If a name is
        provided, it is used (adding '.pkl' if needed); otherwise it is
        presumed this store was created from an existing archive, and
        that archive's name is used.

        In either case, any existing version is backed up before the
        new one is created, with up to two levels of backups (with '%'
        and '%%' suffixes).
        """
        # *** Perhaps this should be called by __del__; or this possibility
        # could be set by a flag on init???

        # If no name given, use the name provided at creation.
        if fname is None:
            if self._fname is None:
                raise ValueError('Need a file name for saving!')
            fname = self._fname
        else:
            if not fname.endswith('.pkl'):
                fname = fname + '.pkl'
        # Gather attributes to store; note new versions supersede old ones.
        attrs = {}
        for name, val in self._pulled.items():
            attrs[name] = val
        # Don't forget to get any archived values not already accessed.
        if self._archive:
            for name in list(self._archive.keys()):
                if name not in self._pulled:
                    attrs[name] = self._archive[name]
        for name, val in self._new.items():
            attrs[name] = val
        # Backup any existing file, and save the attrs to a new file.
        ofile = open(fname, 'wb')
        pkl.dump(attrs, ofile, -1)  # use highest protocol
        ofile.close()

    def contents(self):
        """
        List the names of arrays accessible by this store, including
        both previously archived arrays and newly defined arrays.
        """
        names = []
        if self._archive:
            names.extend(self._archive.keys())
        names.extend(self._new.values())
        return names
