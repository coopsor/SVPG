from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

setup(
    name='svpg',
    version='1.4.0',
    description='A pangenome-based structural variant caller',
    author='henghu',
    author_email='hhengwork@gmail.com',
    license="MIT",
    long_description=readme,
    long_description_content_type='text/markdown',
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        'numpy>=1.26.4',
        'pysam>=0.22',
        'scipy>=1.13.1',
        'pyabpoa>=1.5.4',
        # 'mappy>=2.28',
    ],
    entry_points={
        'console_scripts': [
            'svpg=svpg.main:main',
        ],
    },
    python_requires='>=3.10',
)
