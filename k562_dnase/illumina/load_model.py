import imp
#load the model 
architecture_module=imp.load_source('','bassenji_2019.py')
args=[]
model=architecture_module.getModelGivenModelOptionsAndWeightInits(args)
