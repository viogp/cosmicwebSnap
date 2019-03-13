### e-mails ###

*15th February 2019, from Michaela*
> > I didn't understand the description you sent about the dust attenuation. 
> > Firstable, the tables al*dat: do they have wavelength(nanometers), Al(wavelength) ?  

They are in Angstroem. 

> > If the tables contain that (A(wavelength)=Al and wavelength), why I cannot simply do: F_obs = F_int*10^(-0.4 Al) ? I guess, the Al referes to the extinction curve, which I've usually seen referred as kl. Just to double check with you then, what I have in the tables is an extinction curve normalised  in such a way that what I need to do is: 
> > F_obs = F_int*10^(-0.4 * Av * c2(interpolated at leff)/c2(interpolated at 550nm)) 
> > with c2 being the second column in the file and leff= the wavelenght of the line.  

Exactly, that’s how I would do it!! 

> > In at  Galform, the Av is computed from the properties of the galaxy, through a combination of the properties of its cold gas. In this case, do you think assuming Av=1 is ok? I could play with the fixed value, but my concern here is more the fact that we are fixing Av. I had a look to De Lucia+04, but I didn't find what values of Av are assumed. Could you let me know what do you usually assume?  

I would say that Av=1 is completely fine. That’s what we normally assume in GAEA as well…

*7th August 2018, from Michaela*
I finally generated a SAM-based galaxy catalogue with nebular emission lines from ionised gas due to young stars at z=0.99. The catalogue is based on the FIRE model from my 2016 paper, using, as a start, 1/5 of the Millennium volume (i.e., 100 tree-files out of 500). 

The coupling procedure between SAMs and nebular emission line models is as described in my 2017 paper, but neglecting any contribution from AGN and post-AGB stars. One difference in the coupling procedure, with respect to my paper, is that the ionisation parameter logU is not explicitly calculated from gas density and SFR, but from an empirical relation between logU and the gas metallicity (Fig. 2 of Carton et al. 2017).

You can download the catalogue on this link:
https://www.dropbox.com/s/96211l82lk6qijs/NebCat-H17_FIRE-H16_z10.dat.gz

The format is the following:
xpos [kpc], ypos [kpc], zpos [kpc], log(Mhosthalo) [Msolar], type [0: central, 1:satellite, 2: orphant satellite], log(Mhalo) [Msolar], 
          log(Mstellar) [Msolar], SFR [Msolar/yr], log(Mcold) [Msolar], log(Zcold) [Zsolar], log[C/O], 
          u_mag, g_mag, r_mag, i_mag, z_mag, log(lum_Hb) [Lsolar], log(lum_OIII) [Lsolar], log(lum_Ha) [Lsolar], 
	log(lum_OI) [Lsolar], log(lum_SII) [Lsolar], log(lum_NII) [Lsolar], log(lum_OII) [Lsolar]

0 xpos [kpc]
1 ypos [kpc]
2 zpos [kpc]
3 log(Mhosthalo) [Msolar]
4 type [0: central, 1:satellite, 2: orphan satellite]
5 log(Mhalo) [Msolar] 
6 log(Mstellar) [Msolar]
7 SFR [Msolar/yr]
8 log(Mcold) [Msolar]
9 log(Zcold) [Zsolar]
10 log[C/O] 
11 u_mag  (no dust)
12 g_mag
13 r_mag
14 i_mag
15 z_mag
16 log(lum_Hb) [Lsolar]
17 log(lum_OIII) [Lsolar]
18 log(lum_Ha) [Lsolar]
19 log(lum_OI) [Lsolar]
20 log(lum_SII) [Lsolar]
21 log(lum_NII) [Lsolar]
22 log(lum_OII) [Lsolar]

*30th April 2018, to Michaela*

 I'm looking into 2 aspects:
* how does the treatment of satellites affect the clustering of ELGs? Here is where I'd find interesting to have your GAEA model using either a FIRE like stellar feedback or the PREH that you describe in your paper.
* do ELGs live in filaments?

I'm working with a WMAP7 simulation and I'm focusing only on z=0.83 and z=0.99.

The minimal information I'd need is:
x,y,z, L[OII] - preferibly both attenuated by dust and intrinsic, type (satellite/central), sSFR, stellar mass, mass of the host halo.

It would be good to have also AB magnitudes:
R(Deimos), i(CFHT), g,r,z (decam)

If you don't have the Deimos and CFHT filters, I think SDSS, would be good enough.

I'll need the tables until the paper is published. Afterwards, it might be good to make them publicly available somehow.
Michaela Hirschman 13/1/2019

Regarding your paper, I would love to still contribute to the project, and my last teaching will be done in two weeks so that I am also more free for science again! Let me go through your questions below:
>
>
>
>     * L[OII] is this the doublet or a single line? I'm actually using the OII3727 line.

Yes, I have been giving you the same line OII 3727!
>
>
>     * Could you send me the Lsolar that you assumed, so I can express things in the same units?

As Lsolar I considered 0.02.
>
>
>     I send attached a comparison of the passive fractions at z=0. I hope I can use your model, because it seems to nail the observed passive fraction! 
This is don’t understand, why should it be a problem, when Gaea is matching the observed passive fraction?

>     The problem at the moment is that I do the selections on L[OII] for only the attenuated OII3727 line, so I'll need to access this information. 
So I have the dust attenuation curve from calzetti et al. I attach you below the files, and a pages-file with the description of how to use them. Please let me know if you cannot open the pages-file, then I’ll convert it into a pdf..

Let me know if you have any questions, we can also Skype this week, if you want, just let me know which day/time works best for you.

Hope this helps!

Cheers,
Michaela

al_avcalz == Calzetti Rv = 4.05

al_av31 == MW Rv=3.1