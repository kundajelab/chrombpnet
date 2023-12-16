from setuptools import setup,find_packages


#generate install_requires from requirements.txt file
install_requires=open('requirements.txt','r').read().strip().split('\n')
print(f"install_requires:{install_requires}")


config = {
    'name': 'chrombpnet',
    'author_email': 'anusri @ stanford.edu',
    'license': 'MIT',
    'include_package_data': True,
    'description': 'chrombpnet predicts chromatin accessibility from sequence',
    'download_url': 'https://github.com/kundajelab/chrombpnet',
    'version': '0.1.7',
    'packages': find_packages(),
    'python_requires': '>=3.8',
    'install_requires': install_requires,
    'zip_safe': False,
    'scripts':[
               'chrombpnet/training/models/bpnet_model.py',
               'chrombpnet/training/models/chrombpnet_with_bias_model.py'
    ],
    'entry_points': {'console_scripts': [
        'chrombpnet = chrombpnet.CHROMBPNET:main',
        'print_meme_motif_file = chrombpnet.data.__init__:print_meme_motif_file']}
}

if __name__== '__main__':
    setup(**config)
