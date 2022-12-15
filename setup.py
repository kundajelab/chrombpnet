from setuptools import setup,find_packages


#generate install_requires from requirements.txt file
install_requires=open('requirements.txt','r').read().strip().split('\n')
print(f"install_requires:{install_requires}")


config = {
    'name': 'chrombpnet',
    'author': 'Kundaje lab',
    'author_email': 'anusri @ stanford.edu',
    'license': 'MIT',
    'include_package_data': True,
    'description': 'chrombpnet predicts chromatin accessibility from sequence',
    'download_url': 'https://github.com/kundajelab/chrombpnet',
    'version': '1.3',
    'packages': find_packages(),
    'python_requires': '>=3.8',
    'install_requires': install_requires,
    'zip_safe': False,
    'scripts':['chrombpnet/helpers/make_gc_matched_negatives/make_gc_matched_negatives.sh',
               'chrombpnet/training/models/bpnet_model.py',
               'chrombpnet/training/models/chrombpnet_with_bias_model.py',
               'workflows/tutorial/step1_download_bams_and_peaks.sh',
               'workflows/tutorial/step2_make_bigwigs_from_bams.sh',
               'workflows/tutorial/step3_get_background_regions.sh',
               'workflows/tutorial/step4_train_bias_model.sh',
               'workflows/tutorial/step5_interpret_bias_model.sh',
               'workflows/tutorial/step6_train_chrombpnet_model.sh',
               'workflows/tutorial/step7_interpret_chrombpnet_model.sh',
               'workflows/train_bias_model.sh',
               'workflows/train_chrombpnet_model.sh'
    ],
    'entry_points': {'console_scripts': [
        'chrombpnet_genomewide_gc = chrombpnet.helpers.make_gc_matched_negatives.get_genomewide_gc_buckets.get_genomewide_gc_bins:main',
        'chrombpnet_gc_content_foreground = chrombpnet.helpers.make_gc_matched_negatives.get_gc_content:main',
        'chrombpnet_gc_matched_negatives = chrombpnet.helpers.make_gc_matched_negatives.get_gc_matched_negatives:main',
        'chrombpnet_make_splits = chrombpnet.helpers.make_chr_splits.splits:main',
        'chrombpnet = chrombpnet.CHROMBPNET:main',
        'print_meme_motif_file = chrombpnet.data.__init__:print_meme_motif_file']}
}

if __name__== '__main__':
    setup(**config)
