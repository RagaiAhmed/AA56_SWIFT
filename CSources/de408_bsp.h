/* Constants for the JPL DE421 Chebyshev Ephemeris
 */

#define DE102F 0
#define DE200F 0
#define DE245F 0
#define DE400F 0
#define DE421BSP 0
#define DE408BSP 1
#define BSP_FORMAT 1

double JDb = 1.5; /* starting JED -10020 May 26 */
double JDe = 5376912.5; /* ending JED 10009 May 21 */
double JDi = 16.0; /* days per record varies */
long recsiz = 16; /* bytes per record varies */
long rec0 = 4096L; /* offset to first record varies */

struct ephloc
{
  long first_item_ordinal;
  long last_item_ordinal;
  int days_per_record;
  int items_per_record;
  int total_records;
};

struct ephloc objs[15] = {
  {513,    20116964, 16, 44, 457192},      /* Mercury */
  {20116965, 24460292, 64, 38, 114298},  /* Venus */
  {24460293, 31089580, 32, 29, 228596},  /* EMB */
  {31089581, 34747120, 64, 32, 114298},  /* Mars */
  {34747121, 37033084, 64, 20, 114298},
  {37033085, 39319048, 64, 20, 114298},
  {39319049, 41605012, 64, 20, 114298},
  {41605013, 43890976, 64, 20, 114298},
  {43890977, 46176940, 64, 20, 114298},
  {46176941, 50520268, 64, 38, 114298},  /* Sun */
  {50520269, 88010016, 8, 41, 914384}, /* Moon re EMB */
  {88010017, 125499764, 8, 41, 914384}, /* Earth re EMB */
  {125499765, 125499776, 7315072, 8, 1},  /* Mercury point */
  {125499777, 125499788, 7315072, 8, 1},  /* Venus point */
  {125499789, 125499800, 7315072, 8, 1}  /* Mars point */
#if 0
  {  3, 14, 4},
  {171, 10, 2},
  {231, 13, 2},
  {309, 11, 1},
  {342, 8, 1},
  {366, 7, 1},
  {387, 6, 1},
  {405, 6, 1},
  {423, 6, 1},
  {441, 13, 8},
  {753, 11, 2},
  {819, 10, 4},
  {899, 10, 4}
#endif
};

double au = 1.49597870691e8; /* Kilometers per au */
double emrat = 81.30056; /* Earth/Moon mass ratio */
double radearth = 6378.137;
double clight = 2.99792458e5;
char *dename = "de408.bsp";
