from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

ext_modules = [
    Extension(f"licofage.{pyx}", [f"src/pyx/{pyx}.pyx"])
    for pyx in "addcut dfamin arraysorter".split()
]


class LazyImportBuildExtCmd(build_ext):
    def run(self):
        import numpy as np

        self.include_dirs.append(np.get_include())
        super().run()


setup(
    cmdclass={"build_ext": LazyImportBuildExtCmd},
    ext_modules=ext_modules,
)
