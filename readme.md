# respknt

Reflection matrix approach to computing the seismic response of a cylindrically symmetric medium.

This program was written by George Randall and is based on Kennett's (1983). 

Introduction: http://eqseis.geosc.psu.edu/~cammon/HTML/RftnDocs/resp01.html

Forked from George Randall's code

Reference: 

George E. Randall; Efficient calculation of complete differential seismograms for laterally homogeneous earth models, Geophysical Journal International, Volume 118, Issue 1, 1 July 1994, Pages 245–254, https://doi.org/10.1111/j.1365-246X.1994.tb04687.x

```

1_src: only pick up the useful files.
       no other dependences
2_src: put the subs in src
       can put together/ no need the link file
3_src: add interface.f for python test
       By f2py.
4_rf_respknt_theo:
       compare result of respknt with theo
5_rf_theo_respknt_raysum
       compare result of respknt with theo and raysumy
```