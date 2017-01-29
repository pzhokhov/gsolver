#ifndef _PULSE_SHAPES
#define _PULSE_SHAPES


const char SHAPE_GG[] = "gg";
const char SHAPE_SS[] = "ss";
const char SHAPE_TWOHUMPSGG[]      = "twohumpsgg";
const char SHAPE_TWOBEAMSGG[]      = "twobeamsgg";
const char SHAPE_TBPUMPPROBEGG[]   = "tbpumpprobegg";
const char SHAPE_PUMPPROBEXGG[]    = "pumpprobexgg";
const char SHAPE_BEAMSGG[]         = "beamsgg";
const char SHAPE_BEAMSPROBEGG[]    = "beamsprobegg";
const char SHAPE_BEAMSPROBEEXTGG[] = "beamsprobeextgg";
const char SHAPE_BEAMSPROBEVGG[]   = "beamsprobevgg";
const char SHAPE_BEAMS_DELAYEDPROBE_GG[] = "beams_delayedprobe_gg";
const char SHAPE_BEAMS_DELAYEDPROBE_RANDPHASE[] = "beams_delayedprobegg_randphase";
const char SHAPE_CUSTOMSPECTRUMG[] = "customspectrumg";

void create_gg					(FILE* fid);
void create_beamsgg				(FILE* fid);
void create_beamsprobegg			(FILE* fid);
void create_beamsprobeextgg			(FILE* fid);
void create_beamsprobevgg			(FILE* fid);
void create_beams_delayedprobegg		(FILE* fid);
void create_beams_delayedprobegg_randphase	(FILE* fid);
void create_customspectrumg			    (FILE* fid);



#endif
