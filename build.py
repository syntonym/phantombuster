import numpy as np
import pyarrow as pa

import os

from setuptools.extension import Extension
from Cython.Build import cythonize

def build(setup_kwargs):
    """
    This is a callback for poetry used to hook in our extensions.
    """
    ext_modules = cythonize("src/phantombuster/merge_cython.pyx")

    import pyarrow
    pyarrow.create_library_symlinks()

    for ext in ext_modules:
        # The Numpy C headers are currently required
        ext.include_dirs.append(np.get_include())
        ext.include_dirs.append(pa.get_include())
        ext.libraries.extend(pa.get_libraries())
        ext.library_dirs.extend(pa.get_library_dirs())

        if os.name == 'posix':
            ext.extra_compile_args.append('-std=c++17')
        ext.runtime_library_dirs.append("$ORIGIN/../pyarrow")

        # Try uncommenting the following line on Linux
        # if you get weird linker errors or runtime crashes
        ext.define_macros.append(("_GLIBCXX_USE_CXX11_ABI", "0"))


    setup_kwargs.update(
        {
            # declare the extension so that setuptools will compile it
            "ext_modules": ext_modules,
        }
    )
