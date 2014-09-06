// Read all of the raw UCAC4 data, isolate stars brighter than mvS mv and save
// these to files each containing about a square degree. Save those brighter
// than mvC mv to a candidate list of possible double star primaries.

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// The UCAC4 goes down to 16mv. The fainter you set these limits to, the more
// stars in your final candidate list and square degree regions. Note that
// magnitudes are expressed in millimagnitudes.  12,000 = 12.0mv.
int mvC = 11000, // The minimum brightness of a candidate star.
    mvS = 12000; // The minimum brightness of a star star to save.

// 30" in radians.
const double margin = 3.14159265358979323846 / (180 * 60 * 2);
const double pi = 3.14159265358979323846;

FILE *CAN,    // The stars that might be a new pair's primary.
     *MGN,    // A zone near the margin of the current zone.
     *SQUARE, // The square degree files that stars to be studied are stored in.
     *RAW;    // The UCAC4 raw data.

int cCt = 0,    // The number of candidate stars found  
    starCt = 0; // The number of stars brighter than 14mv that are studied.

// Store this data from a given UCAC4 entry.
typedef struct Candidate_Data {
  char deg[24];// The square degree a candidate is located in.
  double ra;   // Right ascension in radians.
  double dec;  // Declination in radians.
  int dFlg;    // UCAC4 double flag.
  int id;      // UCAC4 zone id.
  int mv;      // Visual magnitude.
  int mvs;     // Source: 0 = APASS, 1 = UCAC4 model. 
  int pmRa;    // Proper motion in right ascension in mas/year.
  int pmDec;   // Proper motion in declination in mas/year.
  int zone;    // UCAC4 zone.
} cData;


// Store this data from a given UCAC4 entry.
typedef struct UCAC4_Data {
  double ra;   // Right ascension in radians.
  double dec;  // Declination in radians.
  int dFlg;    // UCAC4 double flag.
  int id;      // UCAC4 zone id.
  int mv;      // Visual magnitude.
  int mvs;     // Source: 0 = APASS, 1 = UCAC4 model. 
  int pmRa;    // Proper motion in right ascension in mas/year.
  int pmDec;   // Proper motion in declination in mas/year.
  int zone;    // UCAC4 zone.
} uData;

int main(int argc, char** argv) {
  void openRawData(int i), // Close current raw data file and open a new one.
       processRawData();   // Read the raw file and convert it to NA format.
  
  time_t start = time(0);

  // This is my holding directory on my machine.  Adjust it to fit your own.
  // The directory must be empty before we begin.
  system("/bin/rm /science/tmp/*");

  // The primary star candidates are kept here.
  CAN = fopen("/science/tmp/candidates", "a");
  if (CAN == 0) {
    printf("File candidates was not opened!\n");
    exit(0);
  }

  // Parse through all of the stars in the UCAC4.
  for (int i = 1; i < 901; i++) {
    openRawData(i);
    processRawData(i);
  }

  // We're done. Close the still open files.
  fclose(CAN);
  fclose(RAW);
  if (SQUARE) { fclose(SQUARE); }

  time_t end = time(0);
  int delta = (int) (end - start);
  int hr = delta / 3600;
  int min = (delta - (hr * 3600)) / 60;
  int sec = delta % 60;
  printf("Done. Found %d stars. The run took %d:%d:%d.\n",
         starCt, hr, min, sec);

  printf("Found %d candidate stars.\n",cCt); //TEST
}

// Close current raw data file and open a new one.
void openRawData(int i) {

  char df[64],
       fileName[8];

  // Where the raw WDS file live on my machine. Adjust it to point to your own
  // repository.
  static char* rawDir = "/science/astro/data/ucac4/data/";
  df[0] = 0;
  if (RAW) { fclose(RAW); }
 
  sprintf (fileName, "z%03d", i);
  strcat(df, rawDir);
  strcat(df, fileName);

  RAW = fopen(df, "r");
  if (RAW == NULL) {
    printf("Couldn't open %s because\n  %s.\n", df, strerror(errno));
    exit(1);
  }
}

// Read the raw file and convert it to an ASCII format.
// RA and dec are converted to in radians.
void processRawData(int zone) {
  char str[24];     // Generic string.
  char sqDeg[24];   // Name of a square degree zone file.
  int curDec = -99999,
      numLines = 1,
      curRa = -99999,
      sCt = 1;

  unsigned char d[78];

  while (numLines == 1) {
    numLines = fread(&d, 78, 1, RAW);
    if (numLines != 1) { break; }

    int mv = (d[49] * 256) + d[48]; // APASS mv.
    int mvSource = 0;
    if (mv == 20) {
      mv = (d[9] * 256) + d[8];
      mvSource = 1;
    }

    if (mv > mvS) { continue; } // Stars must be brighter than 14mv.
    starCt++;

    // Convert milliarcseconds to radians.
    double raMas = (double) ((d[3] << 24) + (d[2] << 16) + (d[1] << 8) + d[0]);
    double raRad = (double) raMas * pi / (3600000 * 180);

    double decMas = (double) ((d[7] << 24) +(d[6] << 16) + (d[5] << 8) + d[4]);
    double decRad = (((double) decMas / 3600000) - 90) * pi / 180;

    // These "coordinates" identify the file to store the star in.
    // They are in units of integer degrees.
    double raDeg = raRad * cos(decRad) * 180 / pi;
    int ra = (int) raDeg; 
    double decDeg = decRad * 180 / pi;
    int dec = (int) decDeg; 

    uData uStar;
    uStar.ra = raRad;
    uStar.dec = decRad;
    uStar.mv = mv;
    uStar.mvs = mvSource;
    uStar.pmRa = (d[25] * 256) + d[24];
    if (uStar.pmRa > 32767) { uStar.pmRa -= 65536; }
    uStar.pmDec = (d[27] * 256) + d[26];
    if (uStar.pmDec > 32767) { uStar.pmDec -= 65536; }
    uStar.dFlg = d[14];
    uStar.zone = zone;
    uStar.id = sCt;

    if ((curRa != ra) || (curDec != dec)) {
      if (SQUARE) { fclose(SQUARE); }
      sprintf(sqDeg, "/science/tmp/f%d_s%d", ra, (dec + 89));
      SQUARE = fopen(sqDeg, "a");
      curRa = ra;
      curDec = dec;
    }

    fwrite(&uStar, sizeof(uStar), 1, SQUARE);

    if (mv < mvC) {
      // This is a candidate star.
      cData cStar;
      strcpy(cStar.deg, sqDeg);
      cStar.ra = raRad;
      cStar.dec = decRad;
      cStar.mv = mv;
      cStar.mvs = mvSource;
      cStar.pmRa = uStar.pmRa;
      cStar.pmDec = uStar.pmDec;
      cStar.dFlg = d[14];
      cStar.zone = zone;
      cStar.id = sCt;
      fwrite(&cStar, sizeof(cStar), 1, CAN);
      cCt++;
    }

    // About three percent of the stars will be so close to the edge of a
    // file's boundary, they will need to be in the other file as well.
    double delta = fabs(decDeg - (double) dec);
    if (delta < margin) { // Check the southern boundary.
      double sDec = decDeg - 1;
      if (sDec > -90) {
        sDec *= pi / 180;
        // These "coordinates" again identify the file to store the star in.
        int r = (int) (raDeg * cos(sDec));
        int d = (int) (sDec * 180 / pi) + 89;
        sprintf(str, "/science/tmp/f%d_s%d", r, d);
        MGN = fopen(str, "a");
        fwrite(&uStar, sizeof(uStar), 1, MGN);
        fclose(MGN);
      }
    } else if (delta > (1 - margin)) { // Check the northern boundary.
      double sDec = decRad + 1;
      if (sDec < 90) {
        sDec *= pi / 180;
        // These "coordinates" again identify the file to store the star in.
        int r = (int) (raDeg * cos(sDec));
        int d = (int) (sDec * 180 / pi) + 89;
        sprintf(str, "/science/tmp/f%d_s%d", r, d);
        MGN = fopen(str, "a");
        fwrite(&uStar, sizeof(uStar), 1, MGN);
        fclose(MGN);
      }
    }

    delta = fabs(raDeg - (double) ra);
    if (delta < margin) {               // Check the western boundary.
      double sRa = raDeg - 1;
      if (sRa > 0) {
        // These "coordinates" again identify the file to store the star in.
        int r = (int) (sRa * cos(decRad));
        sprintf(str, "/science/tmp/f%d_s%d", r, (dec + 89));
        MGN = fopen(str, "a");
        fwrite(&uStar, sizeof(uStar), 1, MGN);
        fclose(MGN);
      } else {
        sRa = sRa + 360;
        // These "coordinates" again identify the file to store the star in.
        int r = (int) (sRa * cos(decRad));
        sprintf(str, "/science/tmp/f%d_s%d", r, (dec + 89));
        MGN = fopen(str, "a");
        fwrite(&uStar, sizeof(uStar), 1, MGN);
        fclose(MGN);
      }
    } else if (delta > (1 - margin)) {  // Check the eastern boundary.
      double sRa = raDeg + 1;
      if (sRa < 360) {
        sRa *= pi / 180;
        // These "coordinates" again identify the file to store the star in.
        int r = (int) (sRa * cos(decRad));
        sprintf(str, "/science/tmp/f%d_s%d", r, (dec + 89));
        MGN = fopen(str, "a");
        fwrite(&uStar, sizeof(uStar), 1, MGN);
        fclose(MGN);
      } else {
        sRa = (sRa - 360) * (pi / 180);
        // These "coordinates" again identify the file to store the star in.
        int r = (int) (sRa * cos(decRad));
        sprintf(str, "/science/tmp/f%d_s%d", r, (dec + 89));
        MGN = fopen(str, "a");
        fwrite(&uStar, sizeof(uStar), 1, MGN);
        fclose(MGN);
      }
    }
    sCt++;
  }
}
