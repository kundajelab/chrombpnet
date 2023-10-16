import tensorflow as tf

def get_strategy(args):
    # get tf strategy to either run on single, or multiple GPUs
    if vars(args).get('multiGPU') and args.multiGPU:
        # run the model in "data parallel" mode on multiple GPU devices (on one machine).
        strategy = tf.distribute.MirroredStrategy()
        print('Number of GPU devices: {}'.format(strategy.num_replicas_in_sync))
    else:
        strategy = tf.distribute.get_strategy()
        print('Single GPU device')

    return strategy
