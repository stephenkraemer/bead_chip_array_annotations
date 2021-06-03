from setuptools import find_packages, setup

setup(
    name='mouse_methylation_bead_chip',
    version='0.1',
    author='Stephen Kraemer',
    author_email='stephenkraemer@gmail.com',
    license='MIT',
    packages = ['mouse_methylation_bead_chip'],
    package_data={
        '': ['*.R', '*.snakefile', '*.yml', '*.yaml', '*.sh', '*.smk', '*.rules'],
    },
)
