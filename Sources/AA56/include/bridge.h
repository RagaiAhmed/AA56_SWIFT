#ifndef AA56_SWIFT
#define AA56_SWIFT

// Struct to store polar values
struct Polar {
  double lon;  // Longitude in radians
  double dec;  // Declination in radians
  double r;  // Distance in AU
};

// Initializes the calculator with location on earth
/*
    path: aa.ini file path to init the calculator
        notice to use
    returns 0 if init successfully

    Note: It is prefered to be called before each calc
*/
#ifdef __cplusplus
extern "C"
{
#endif
int initCalc(char path[]);
#ifdef __cplusplus
}
#endif


// Calculates the polar coordinates of a given planet
/*
JulianDate: time in days (Julian Date)

Planet:
0 -> Sun
1 -> Mercury
2 -> Venus
3 -> Moon
4 -> Mars
5 -> Jupiter
6 -> Saturn
7 -> Uranus
8 -> Neptune
9 -> Pluto
*/
#ifdef __cplusplus
extern "C"
{
#endif
struct Polar calcPolar(double JulianDate, int planet);
#ifdef __cplusplus
}
#endif
// Calculates the polar coordinates of a given orbit/star
/*
JulianDate: time in days (Julian Date)

index:
88 -> Orbit
99 -> Star

path: File path of the orbit or star (refer to aa56 readme for more details)
orbit.cat & star.cat are examples
*/
#ifdef __cplusplus
extern "C"
{
#endif
struct Polar calcPolarPath(double JulianDate, int index, char path[]);
#ifdef __cplusplus
}
#endif
#endif // AA56_SWIFT
