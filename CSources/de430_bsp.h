/* Constants for the JPL DE430 Chebyshev Ephemeris
 */

#define DE102F 0
#define DE200F 0
#define DE245F 0
#define DE400F 0
#define DE430BSP 1
#define BSP_FORMAT 1

double JDb = 2288184.5; /* starting JED J1550.055 */
double JDe = 2689976.5; /* ending JED  J2650.0520 */
double JDi = 16.0; /* days per record varies */
long recsiz = 16; /* bytes per record varies */
long rec0 = 5120; /* 4096L; */ /* offset to first record varies */

struct ephloc
{
  long first_item_ordinal;
  long last_item_ordinal;
  int days_per_record;
  int items_per_record;
  int total_records;
};

struct ephloc objs[15] = {
  {    641, 2210500,  8, 44, 50224},  /* Mercury */
  {2210501, 3014088, 16, 32, 25112},  /* Venus */
  {3014089, 4043684, 16, 41, 25112},  /* EMB */
  {4043685, 4483148, 32, 35, 12556},  /* Mars */
  {4483149, 4809608, 32, 26, 12556},  /* Jupiter */
  {4809609, 5098400, 32, 23, 12556},  /* Saturn */
  {5098401, 5349524, 32, 20, 12556},  /* Uranus */
  {5349525, 5600648, 32, 20, 12556},  /* Neptune */
  {5600649, 5851772, 32, 20, 12556},  /* Pluto */
  {5851773, 6730696, 16, 35, 25112},  /* Sun */
  {6730697, 10849068, 4, 41, 100448}, /* Moon re EMB */
  {10849069,14967440, 4, 41, 100448}, /* Earth re EMB */
  {14967441,14967452, 401792, 8, 1},  /* Mercury point */
  {14967453,14967464, 401792, 8, 1},  /* Venus point */
};

double au = 149597870.700; /* Kilometers per au */
double emrat =  0.813005690741906200e02; /* Earth/Moon mass ratio */
double radearth = 6378.1363;
double clight = 2.99792458e5;
char *dename = "de430.bsp";
