import numpy as np
import pickle


def save_df(config):
    if hasattr(config, 'main_df'):
        np.savez_compressed(config.config_path + 'main_df', data=config.main_df)
        np.savez_compressed(config.config_path + 'main_df_classes', data=config.main_df_classes)
        pickle_out = open(config.config_path + 'gene_col_dict.pkl', "wb")
        pickle.dump(config.gene_col_dict, pickle_out)
        pickle_out.close()
