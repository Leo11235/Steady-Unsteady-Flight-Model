# ensure that all the PROPEP functions return the right values
from pyPROPEP_simple import call_PROPEP, pyPROPEP_interpolation_lookup, get_chamber_properties_with_partials


class PROPEP_test:
    def __init__(self, OF, CP):
        self.OF = OF
        self.CP = CP
        self.compare_different_functions()
    
    def compare_different_functions(self):
        # check if call_PROPEP, pyPROPEP_interpolation_lookup, get_chamber_properties_with_partials's outputs all match to a reasonable degree
        base_chamber_temp, base_molar_weight, base_heat_ratio = call_PROPEP(self.OF, self.CP)
        interp_x = pyPROPEP_interpolation_lookup(self.OF, self.CP)
        (interp_chamber_temp, interp_molar_weight, interp_heat_ratio) = interp_x
        partial_chamber_temp, partial_molar_weight, partial_heat_ratio, dT_dOF, dW_dOF, dT_dp, dW_dp = get_chamber_properties_with_partials(self.OF, self.CP)
        
        print(f'\n==================== For inputs OF = {self.OF} and CP = {self.CP} ====================')
        print("Base PROPEP            Interpolation                Partials")
        print(base_chamber_temp, interp_chamber_temp, partial_chamber_temp)
        print(base_molar_weight, interp_molar_weight, partial_molar_weight)
        print(base_heat_ratio, interp_heat_ratio, partial_heat_ratio)
        

test0 = PROPEP_test(1, 14.7)
test1 = PROPEP_test(5, 600)
test2 = PROPEP_test(10, 1000)
