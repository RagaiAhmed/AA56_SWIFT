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
int initCalc(char path[]);


// Calculates the polar coordinates of a given planet
/*
UnixTimeStamp: time in seconds

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
struct Polar calcPolar(double UnixTimeStamp, int planet);

// Calculates the polar coordinates of a given orbit/star
/*
UnixTimeStamp: time in seconds

index:
88 -> Orbit
99 -> Star

path: File path of the orbit or star (refer to aa56 readme for more details)
orbit.cat & star.cat are examples
*/
struct Polar calcPolarPath(double UnixTimeStamp, int index, char path[]);

#endif // AA56_SWIFT
