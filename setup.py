from setuptools import setup, find_packages

setup(
    name='svpg',
    version='0.1.0',
    description='A pangenome-based structural variant caller',
    author='henghu',
    author_email='hhengwork@gmail.com',
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        'numpy>=1.26.4',
        'pysam>=0.22',
    ],
    entry_points={
        'console_scripts': [
            'svpg=svpg.main:main',
        ],
    },
    python_requires='>=3.10',
)
