import os
from setuptools import setup, Command

class CleanCommand(Command):
    user_options = []
    def initialize_options(self):
        self.cwd = None
    def finalize_options(self):
        self.cwd = os.getcwd()
    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd
        os.system ('rm -rf ./build ./dist ./*.pyc ./*.tgz ./*.egg-info')

setup(
    name = "update_gene_identifiers",
    version = "0.99.0",
    author = "Yoann Pradat",
    author_email = "yoann.pradat@centralesupelec.fr",
    install_requires = [
        "numpy",
        "pandas",
        "openpyxl",
    ],
    cmdclass={
        'clean': CleanCommand
    }
)

