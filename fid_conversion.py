import nmrglue as ng
from matplotlib import pyplot as plt
import numpy as np

vdic, vdata = ng.varian.read('/Users/leoparsons/Desktop/Coding_Projects/Python_Projecs/NMR_Processing/LP3-049-21hr/', '/Users/leoparsons/Desktop/Coding_Projects/Python_Projecs/NMR_Processing/LP3-049-21hr/fid', '/Users/leoparsons/Desktop/Coding_Projects/Python_Projecs/NMR_Processing/LP3-049-21hr/procpar')

C = ng.convert.converter()
C.from_varian(vdic, vdata)
pdic, pdata = C.to_pipe()

# #Data processing
ng.proc_base.fft(pdata)
ng.proc_base.mir_left(pdata)
ng.proc_base.neg_left(pdata)
ng.proc_bl.sol_sine(pdata)

# #Data transformation
dic,data = ng.pipe_proc.em(pdic,pdata)
dic,data = ng.pipe_proc.zf(dic,data,auto=True)
dic,data = ng.pipe_proc.ft(dic,data,auto=True)
dic,data = ng.pipe_proc.di(dic,data)
data = ng.process.proc_bl.base(data, nl=list(np.arange(0, 32768, 1, dtype=int)))
data = ng.proc_base.rev(data)
data = ng.proc_autophase.autops(data, "acme")



#Plot data


# Phasing. Use the autophase.manula_ps to mess with the phasing then you can
# input those values into proc_base.ps function below

# p0, p1 = ng.proc_autophase.manual_ps(ffted_nrmglue)
# phased_data = ng.proc_base.ps(ffted_nrmglue, p0=88, p1=-282)
# # exp_data = ng.process.proc_base.ps_exp(ffted_nrmglue, p0=600, tc=4, inv=True)
# real_phased_data = ng.proc_base.di(phased_data)
#
#
unit_conv = ng.fileio.fileiobase.unit_conversion(size=np.size(data), cplx=False, sw=35000, obs=1, car=1)
data = unit_conv.ppm(data)

plt.plot(data)
plt.show()

# integ_data_1 = ng.integration.integrate(real_phased_data, unit_conv, [14000, 14200])
# integ_data_2 = ng.integration.integrate(real_phased_data, unit_conv, [15800, 16090])
# integ_std = ng.integration.integrate(real_phased_data, unit_conv, [22670, 21470])
#
# ref_integ_data_1 = integ_data_1 / integ_std
# ref_integ_data_2 = integ_data_2 / integ_std
#
# print(ref_integ_data_1)
# print(ref_integ_data_2)
#
# print("ok")
# plt.plot(phased_data)
# print("ok")
# plt.show()


