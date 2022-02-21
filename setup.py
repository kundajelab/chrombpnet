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
    'python_requires': '>=3.6',
    'install_requires': install_requires,
    'dependency_inks': dependency_links,
    'zip_safe': False,
    'scripts':['chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh',
               'chrombpnet/helpers/preprocessing/bam_to_bigwig.sh',
               'step1_download_bams_and_peaks.sh',
               'step2_make_bigwigs_from_bams.sh',
               'step3_get_background_regions.sh',
               'step4_train_bias_model.sh',
               'step5_interpret_bias_model.sh',
               'step6_train_chrombpnet_model.sh',
               'step7_interpret_chrombpnet_model.sh'
    ],
    'entry_points': {'console_scripts': [    'chrombpnet_hyperparams = chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams:main',
                                             'chromobpnet_bias_hyperparams = chrombpnet.helpers.hyperparameters.find_bias_hyperparams:main',
                                             'chrombpnet_train = chrombpnet.training.train:main',
                                             'chrombpnet_predict = chrombpnet.training.predict:main',
                                             'chrombpnet_metrics = chrombpnet.training.metrics:main',
                                             'chrombpnet_deepshap = chrombpnet.evaluation.interpret.interpret:main',
                                             'chrombpnet_modisco = chrombpnet.evaluation.modisco.run_modisco:main',
                                             'chrombpnet_marginal_footprints = chrombpnet.evaluation.marginal_footprints.marginal_footprinting:main',
                                             'chrombpnet_pwm_from_bigwig = chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig:main',
                                             'chrombpnet_srcdir = chrombpnet.get_package_dir:main']}}

if __name__== '__main__':
    setup(**config)
