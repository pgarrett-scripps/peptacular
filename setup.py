from setuptools import setup

# Prior to pip v21.1, a setup.py script was required to be compatible with development mode. With late versions of pip,
# projects without setup.py may be installed in this mode.
setup(
    include_package_data=True,
    package_dir={"": "src"},
)
