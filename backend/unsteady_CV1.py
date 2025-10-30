# Tank control volume: computes [n_v, n_l, T_T] for each timestep

def tank_CV1(liquid_oxidizer_volume):
    if liquid_oxidizer_volume > 0: # liquid blowdown
        return lquid_oxidizer_blowdown()
    elif liquid_oxidizer_volume == 0: # gaseous blowdown
        return gaseous_oxidizer_blowdown()
    else:
        print("Error in CV1: liquid_oxidizer_volume neither >0 or =0")

def lquid_oxidizer_blowdown():
    return

def gaseous_oxidizer_blowdown():
    return