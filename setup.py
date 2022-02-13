from setuptools import setup,find_packages

#generate install_requires and dependency_links from requirements.txt file
reqs=open('requirements.txt','r').read().strip().split('\n')
install_requires=[]
dependency_links=[]
for line in reqs:
    if line.startswith('git'):
        dependency_links.append(''.join(line.split('+')[1::]))
    else:
        install_requires.append(line)

print(f"dependency_links:{dependency_links}")
print(f"install_requires:{install_requires}")

config = {
    'name': 'chrombpnet',
    'author': 'Kundaje lab',
    'author_email': 'anusri @ stanford.edu',
    'license': 'MIT',
    'include_package_data': True,
    'description': 'chrombpnet predicts chromatin accessibility from sequence',
    'download_url': 'https://github.com/kundajelab/chrombpnet',
    'version': '0.1',
    'packages': ['chrombpnet'],
    'install_requires': install_requires,
    'dependency_inks': dependency_links,
'entry_points': {'console_scripts': [    'chrombpnet_get_params = chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams:main',
                                         'chrombpnet_train = chrombpnet.training.train:main',
                                         'chrombpnet_predict = chrombpnet.training.predict:main',
                                         'chrombpnet_metrics = chrombpnet.training.metrics:main',
                                         'chrombpnet_deepshap = chrombpnet.evaluation.interpret.interpret:main',
                                         'chrombpnet_modisco = chrombpnet.evaluation.modisco.run_modisco:main',
                                         'chrombpnet_marginal_footprints = chrombpnet.evaluation.marginal_footprints.marginal_footprinting:main']}}

if __name__== '__main__':
    setup(**config)
