// Using the candidate pairs and regions created by mkUCAC4_Regions, look for
// pairs of stars that are :
//   -> Within XXX" of each other.
//   -> Within dMv mv of each other.
//   -> Not within XXX" of a WDS pair.
// An HTML format list of unlisted pairs will be created.

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

const double pi = 3.14159265358979323846;

FILE *CAN,    // The stars that might be a new pair's primary.
     *MGN,    // A zone near the margin of the current zone.
     *SQUARE; // The square degree files that stars are stored in.

// UCAC4 magnitudes are expressed in thousandth of a magnitude.
// Changing these parameters will greatly affect the number of unlisted pairs
// found.
int dMv = 4000,   // The maximum magnitude difference between the two stars.
    maxF = 65536, // The maximum number of stars this program will find.
    maxSep = 30,  // The maximum separation in arc seconds.
    minPM = 5,    // The minimum combined proper motion for a pair.
    minSep = 2,   // The minimum separation in arc seconds.
    mvC = 12000,  // The minimum magnitude for a primary candidate star.
    mvS = 13000,  // The minimum magnitude for a secondary star.
    pmR = 2;      // The minimum proper motion to delta proper motion.

int wdsCt = 0,  //TEST
    wdsOut = 0; //TEST

// Stars within a box of XXX arc seconds centered on the primary candidate will
// be considered as companions of the candidate.
const double XXX = 3.14159265358979323846 / (180 * 60 * 2); // 30" in radians.

// Store this data from a given Candidate entry.
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

// Data from a given UCAC4 entry.
typedef struct UCAC4_Data {
  double ra;  // Right ascension in radians.
  double dec; // Declination in radians.
  int dFlg;   // UCAC4 double flag.
  int id;     // UCAC4 zone id.
  int mv;     // Visual magnitude.
  int mvs;    // Source: 0 = APASS, 1 = UCAC4 model. 
  int pmRa;   // Proper motion in right ascension in mas/year.
  int pmDec;  // Proper motion in declination in mas/year.
  int zone;   // UCAC4 zone.
} uData;

// We only need position from the WDS. No candidate pairs within 30" of a WDS
// pair are considered viable.
typedef struct WDS_Data {
  double ra;   // Right ascension in radians.
  double dec;  // Declination in radians.
} wdsD;

wdsD wds[128880]; // The 128880 WDS stars.
int wdsIndex[64]; // Index into the wds array.

int main(int argc, char** argv) {

  void r2ra(char*, double ra); // Convert radian ra to hms ra.
  void r2dec(char*, double d); // Convert radian dec to dms dec.

  int ckWDS(double east,  // Check to see if there is a WDS pair within 
            double north, // these boundaries.
            double south,
            double west);
  void readWDS(void);     // Read the coordinates of all precise WDS pairs.

  time_t start = time(0);

  CAN = fopen("/science/tmp/candidates", "r");
  if (CAN == 0) {
    printf("File candidates was not opened!\n");
    exit(0);
  }

  FILE* NEW = fopen("/work/glxy/tmp/unlistedPairs.html", "w");
  if (NEW == 0) {
    printf("File unlistedPairs was not opened!\n");
    exit(0);
  }
  fprintf(NEW, "\n<!DOCTYPE html PUBLIC Content-type: text/html>\n<HTML>"
          "<BODY BGCOLOR=navy TEXT=white><CENTER>\n"
          "<TITLE>Non WDS pairs</TITLE>\n<H2>Non WDS pairs.</H2><BR>\n"
          "<TABLE BORDER=8><TR><TD>RA Dec</TD><TD>mv</TD><TD>mv src</TD>"
          "<TD>mvb</TD><TD>mvb src</TD><TD>&rho;\"</TD><TD>Double<BR>Flag</TD>"
          "<TD>Primary<BR>PM in RA</TD><TD>Primary<BR>PM in Dec</TD>"
          "<TD>Secondary<BR>PM in RA</TD><TD>Secondary<BR>PM in Dec</TD>"
          "<TD>A UCAC4 id</TD><TD>B UCAC4 id</TD><TD>Comments</TD></TR>\n");

  readWDS(); // Load in the WDS.

  // Open and read the candidates list.
  char sqDeg[24]; // The current square degree being examined.

  int cCt = 0, // Count the candidates as they're checked.
      lCt = 1, // Number of lines read from the candidate list of stars.  
      pCt = 0; // Number of unlisted pairs found.

  int t = 0; //TEST

  cData cStar; // A candidate star.
  while (lCt == 1) {
    lCt = fread(&cStar, sizeof(cData), 1, CAN);
    if (lCt != 1) { break; }

    if (cStar.mv > mvC) { continue; }

    // Stars within this box are considered candidates for pairs.
    double east = cStar.ra + XXX,
           north = cStar.dec + XXX,
           south = cStar.dec - XXX,
           west = cStar.ra - XXX;

    int ckLines = 1;  // Number of lines read from a square degree of stars.

    // Open this square if it's not already open.
    if (strcmp(cStar.deg, sqDeg) != 0) {
      if (SQUARE) { fclose(SQUARE); }
      SQUARE = fopen(cStar.deg, "r");
      if (SQUARE == NULL) {
        printf("Failed to open square degree %s.\n", cStar.deg);
        exit(0);
      }
      strcpy(sqDeg, cStar.deg);

    } else {
      rewind(SQUARE); // Start searching at the start of the file.
    }

    uData ckSt; // A star to check aganist the candidate star.
    while (ckLines == 1) {
      // Is the star in the box, and is it not the same as the candidate?
      ckLines = fread(&ckSt, sizeof(uData), 1, SQUARE);
      if (ckLines != 1) { break; }

      if (ckSt.mv > mvS) { continue; }

      if ((ckSt.ra > west) && (ckSt.ra < east) &&
          (ckSt.dec < north) && (ckSt.dec > south)) {

        // Don't pick the same star as the candidate!
        // if (cStar.id == ckSt.id) { continue; } RESTORE!
        if (cStar.id == ckSt.id) { continue; t++; } //TEST

        // The primary should outshine the secondary. 
        if (cStar.mv >= ckSt.mv) { continue; }

        // Stars need to be within 1mv of each other.
        if (abs(ckSt.mv - cStar.mv) > dMv) { continue; }

        // The proper motion of the secondary must not be zero.
        if ((ckSt.pmRa == 0) && (ckSt.pmDec == 0))  { continue; }

        // The stars neet to be within maxSep arc seconds of each other, but not
        // within minSep.
        double dr = cStar.ra - ckSt.ra;
        double dd = cStar.dec - ckSt.dec;
        double sep = (sqrt((dr * dr) + (dd * dd))) * 180 * 3600 / pi;
        if ((sep < minSep) || (sep > maxSep)) { continue; }

        // A possible double star.  Check the proper motions.
        double pmR = (cStar.pmRa + ckSt.pmRa) / 2,
               pmD = (cStar.pmDec + ckSt.pmDec) / 2,
               pm = sqrt((pmR * pmR) + (pmD * pmD)); 

        // The proper motion should be more than minPM milliarcseconds/yr.
        if (pm < minPM)  { continue; }

        int rDel = (cStar.pmRa - ckSt.pmRa) / 2;
        int dDel = (cStar.pmDec - ckSt.pmDec) / 2;
        double pmDel = sqrt((rDel * rDel) + (dDel * dDel)); 

        if ((pm / pmDel) > pmR) {
          // This looks like a good candidate.  Is is already in the WDS?
          int w = ckWDS(east, north, south, west);
          if (w) {
            char cStr[8];
            if (! cStar.mvs) { sprintf(cStr, "APASS"); }
            else { sprintf(cStr, "UCAC4_M"); }

            char ckStr[8];
            if (! ckSt.mvs) { sprintf(ckStr, "APASS"); }
            else { sprintf(ckStr, "UCAC4_M"); }
    
            char r[16];
            r2ra(r, cStar.ra);
            char d[16];
            r2dec(d, cStar.dec);

            fprintf(NEW, "<TR><TD>%s %s</TD><TD>%d</TD><TD>%s</TD>" 
            "<TD>%d</TD><TD>%s</TD><TD>%5.2f</TD><TD>%d</TD>" 
            "<TD>%d</TD><TD>%d</TD>" 
            "<TD>%d</TD><TD>%d</TD>" 
            "<TD>%d %d</TD><TD>%d %d</TD>" 
            "<TD><CENTER>-</CENTER></TD></TR>\n",
            r, d, cStar.mv, cStr, ckSt.mv, ckStr, sep, cStar.dFlg,
            cStar.pmRa, cStar.pmDec, ckSt.pmRa, ckSt.pmDec,
            cStar.zone, cStar.id, ckSt.zone, ckSt.id);

            // printf("%d,%d ", cStar.id, ckSt.id); //TEST

            if (pCt > maxF) { ckLines = 0; }
            pCt++;
          }
        }
      }
    }
    cCt++;
  }

  fprintf(NEW, "\n</BODY></HTML>\n");
  fclose(CAN);
  fclose(NEW);
  if (SQUARE) { fclose(SQUARE); }

  time_t end = time(0);
  int delta = (int) (end - start);
  int hr = delta / 3600;
  int min = (delta - (hr * 3600)) / 60;
  int sec = delta % 60;

  // printf("Found %d unlisteds. The run took %d:%d:%d.\n",
  //        (pCt / 2), hr, min, sec);

  printf("Found %d unlisteds. WdsCt: %d. t: %d. The run took %d:%d:%d.\n",
         (pCt / 2), wdsCt, t, hr, min, sec);
}

// Check to see if there is a WDS pair within these boundaries.
int ckWDS(double east,
          double north,
          double south,
          double west) {

  int i = wdsIndex[(int) east * 10] - 10;
  if (i < 0) { i = 0; }


  while (1) {
    if ((wds[i].ra < east) && (wds[i].ra > west) &&
        (wds[i].dec < north) && (wds[i].dec > south)) {
      // There's a WDS pair here, so this one's not unlisted.

      wdsCt++; //TEST

      return 0;
    } else if (wds[i].ra > east) {
      // This one looks good.  Save it!

      wdsOut++; //TEST
     
      return 1;
    }
    i++;
  }
}

// Convert radian dec to d:m:s.s dec.
void r2dec(char* str,
           double dec) {
  char sign[8];
  if (dec < 0) {
    sprintf(sign, "-");
    dec *= -1;
  } else {
    sprintf(sign, "+");
  }

  double deg = dec * 180 / pi; // Convert radians to degrees.
  int d = (int) deg;

  double min = (deg - (double) d) * 60;
  int m = (int) min;

  double sec = (min - (double) m) * 60;
  int s = (int) sec;
  double fraction = sec - (double) s;
  int f = (int) (fraction * 100);

  sprintf(str, "%s%d:%d:%d.%d", sign, d, m, s, f);
}

// Convert radian ra to h:m:s.s ra.
void r2ra(char * str,
          double rra) {
  double hr = rra * 12 / pi;
  int h = (int) hr;

  double min = (hr - (double) h) * 60;
  int m = (int) min;

  double sec = (min - (double) m) * 60;
  int s = (int) sec;
  double fraction = sec - (double) s;
  int f = (int) (fraction * 100);

  sprintf(str, "%d:%d:%d.%d", h, m, s, f);
}

// Read the most recent version of the WDS catalog, and only save the high
// precision coordinates.
void readWDS(void) {

  double str2d(char* str, // Convert the characters in str at
               int loc,   // location loc
               int len);  // len characters long to a double.

  const double d2r = pi / 180; // Converts degrees to radians.
  const double h2r = pi / 12;  // Converts hours to radians.

  // You'll need to change this to the directory where you keep the WDS data.
  FILE* WDS = fopen("/work/glxy/wdsTemp/wdsPrecisionCoordinates", "r");
  if (WDS == 0) {
    printf("The WDS catalog was not opened!\n");
    exit(0);
  }
  char line[24]; // A line of the WDS.
  int lCt = 1;   // The number of lines read by fread.
  int oldRa = 1; // The penultimate value of tenths of RA parsed.
  int wCt = 0;   // The number of high precision pairs found.
  wdsD w;        // The current WDS pair being parsed.
  while (lCt == 1) {
    lCt = fread(&line, 19, 1, WDS);
    if (lCt != 1) { break; }

    // If the coordinates do not exist, go on to the next pair.
    if (line[125] == 32) { continue; }

    // Convert the coordinates to radians and save them.
    double d, h, m, s, f; // Degrees, hours, minutes, seconds.

    h = str2d(line, 0, 2);         // RA Hours.
    m = str2d(line, 2, 2);         // RA Minutes.
    s = str2d(line, 4, 2);         // RA Seconds.
    f = (str2d(line, 7, 2) / 100); // Fractional RA seconds.
    s += f;                        // Add fractional RA seconds.
    w.ra = (h + (m / 60) + (s / 3600)) * h2r;

    d = str2d(line, 10, 2);          // Dec minutes.
    m = str2d(line, 12, 2);          // Dec minutes.
    s = str2d(line, 14, 2);          // Dec seconds.
    s += (str2d(line, 17, 1) / 10);  // Fractional Dec seconds.
    w.dec = (d + (m / 60) + (s / 3600)) * d2r;

    if (line[9] == 45) { w.dec *= -1; } // 45 == "-".

    // Save the coordinates.
    int r = (int) (w.ra * 10);
    if (r > oldRa) {
      oldRa = r;
      wdsIndex[r] = wCt;
    }
    wds[wCt++] = w; 
  }
}

// Convert the characters in str at location loc, len characters long
// to a double.
double str2d(char* str,
             int loc,
             int len) {
  int n = 1,
      res = 0;
  loc--;
  for (int i = (loc + len); i > loc; i--) {
    // Multiply the existing number by 10, convert the next one to an integer
    // from ASCII, and add it to the existing number.
    res += n * (int) (str[i] - 48);
    n *= 10; 
  }
  return (double) res;
}


  // double foundDec[maxF], // RA and Dec are uses as hashes to
  //        foundRa[maxF];  // avoid duplicate listings.
        // Have we already found this pair?
        // int i = pCt - 20;
        // if (i < 0) { i = 0; }
        // for (int j = i; j < pCt; j++) {
        //   if ((ckSt.ra == foundRa[j]) && (ckSt.dec == foundDec[j])) {
        //     continue;
        //   }
        // } 
            // foundRa[pCt] = cStar.ra;
            // foundDec[pCt++] = cStar.dec;
            // foundRa[pCt] = ckSt.ra;
            // foundDec[pCt++] = ckSt.dec;
