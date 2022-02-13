from setuptools import setup,find_packages
import distutils.text_file
from Pathlib import Path
from typing import List

def _parse_requirements(filename: str) -> List[str]:
    """Return requirements from requirements file."""
    return distutils.text_file.TextFile(filename=str(Path(__file__).with_name(filename))).readlines()
            

config = {
    'include_package_data': True,
    'description': 'chrombpnet predicts chromatin accessibility from sequence',
    'download_url': 'https://github.com/kundajelab/chrombpnet',
    'version': '0.1',
    'packages': ['chrombpnet'],
    'setup_requires': [],
    'install_requires': _parse_requirements('requirements.txt')
    'scripts': [],
'entry_points': {'console_scripts': [    'chrombpnet_get_params = chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams:main',
                                         'chrombpnet_train = chrombpnet.training.train:main',
                                         'chrombpnet_predict = chrombpnet.training.predict:main',
                                         'chrombpnet_metrics = chrombpnet.training.metrics:main',
                                         'chrombpnet_deepshap = chrombpnet.evaluation.interpret.interpret:main',
                                         'chrombpnet_modisco = chrombpnet.evaluation.modisco.run_modisco:main',
                                         'chrombpnet_marginal_footprints = chrombpnet.evaluation.marginal_footprints.marginal_footprinting:main']},
    'name': 'chrombpnet'
}

if __name__== '__main__':
    setup(**config)
