import SFLR_TF_models

def create_ROM_model(bode_opts, TMM_model, Gth, Ga):
    ROM_ATFB_model = SFLR_TF_models.G_th_G_a_TF(bode_opts, \
                                                [], \
                                                TMM_model, \
                                                Gth=Gth, \
                                                Ga=Ga, \
                                                ffit=None, \
                                                label='ROM')
    ROM_ATFB_model.load_params('decent_ROM_params_09_09_10.pkl')
    return ROM_ATFB_model
