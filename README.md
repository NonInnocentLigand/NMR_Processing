# NMRProcessor

A package to extend the functionality of the [NMRglue](https://www.nmrglue.com/) package.

---

## Curent Functionality

SingleNMRInitializer

- process, integrate and plot 1d NMR spectra from fid and procpar files

MultipleNMRInitializer

- process integrate and plot multiple 1d NMR spectra
- plot time resolved species concentration using an internal standard signal and a reference signal from the species of intrist

---

### Future Functionality

- fit plots of time resolved with multivariable regression

- highlight integrated peaks on plots

---

### Example Outputs

3d plots of multiple 1d NMR spectra
<p>
    <img src = '1-1-0.2 - 6-1-377 to_3d_plot() fig.png' width='250' />
    <img src = '1-2-0.15 - 5-2-180 to_3d_plot() fig.png' width='250' />
</p>

plots of concentration vs time
<p>
    <img src = '1-1-0.2 to 6-1-377 integrations at ppm = 1.42.png' width='250' />
    <img src = '1-2-0.15 to 5-2-180 integrations at ppm = 1.25.png' width='250' />
</p>