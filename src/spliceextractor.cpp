#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <vector>
#include <cstring>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <math.h>
#include "htslib/sam.h"

#include "arg_parse.h"

#define CMATCH  0
#define CINS    1
#define CDEL    2
#define CREF_SKIP   3
#define CSOFT_CLIP  4
#define CHARD_CLIP  5
#define CPAD    6
#define CEQUAL  7
#define CDIFF   8

// parse cigar and return whether the read is spliced or not. If spliced the coords vector is populated with one or more of the following structure:
// 1. left overhang
// 2. intron length
// 3. right overhang
int process_cigar(bam1_t *al,std::vector<int> &coords){
    bool spliced = false;
    bool first = true;
    int lol=0,rol=0;
    for (uint8_t c=0;c<al->core.n_cigar;++c){
        uint32_t *cigar_full=bam_get_cigar(al);
        int opcode=bam_cigar_op(cigar_full[c]);
        int length=bam_cigar_oplen(cigar_full[c]);

        if (opcode==BAM_CINS){
            return false;
        }
        if (opcode==BAM_CDEL){
            return false;
        }
        if (opcode==BAM_CSOFT_CLIP){
            continue;
        }
        if (opcode == BAM_CMATCH) {
            rol+=length;
            lol+=length;
        }
        if (opcode == BAM_CREF_SKIP){
            spliced = true;
            if(!first){// if a previous entry existed - write its right overhang
                coords.push_back(rol);
            }
            coords.push_back(lol);
            coords.push_back(length);
            lol = 0;
            rol = 0;
            first = false;
        }
    }
    if(spliced){
        coords.push_back(rol);
    }
    return spliced;
}

void add_num_splice_tag(bam1_t *al, int numN){
    uint8_t* ptr_nh_1=bam_aux_get(al,"ZZ");
    if(ptr_nh_1){
        bam_aux_del(al,ptr_nh_1);
    }
    bam_aux_append(al,"ZZ",'i',4,(uint8_t*)&numN);
}

enum Opt {INPUT     = 'i',
          OUTPUT    = 'o',
          NUM_SPLICE= 'n',
          OVERHANG  = 'v'};

// spliceextractor in.bam/cram/sam out_basename max_overhang min_num_introns

int main(int argc, char** argv) {

    ArgParse args("spoon help");
    args.add_string(Opt::INPUT,"input","","input alignment in BAM/CRAM format",true);
    args.add_string(Opt::OUTPUT,"output","","basename of the output files",true);
    args.add_int(Opt::NUM_SPLICE,"num-splice",2,"minimum number of splice junctions to extract reads into a separate BAM file",true);
    args.add_int(Opt::OVERHANG,"overhang",5,"maximum overhang to consider for output",true);

    if(argc == 1){
        std::cerr<<args.get_help()<<std::endl;
        exit(-1);
    }

    if(strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(-1);
    }

    args.parse_args(argc,argv);

    std::cout<<"working with: "<<args.get_string(Opt::INPUT)<<std::endl;
    std::cout<<"writing output to: "<<args.get_string(Opt::OUTPUT)<<std::endl;

    std::string out_base(args.get_string(Opt::OUTPUT));
    std::string out_agg_fname = out_base, out_sam_fname = out_base, out_sam_manyN_fname = out_base;
    out_agg_fname.append(".agg");
    out_sam_fname.append(".bam");
    out_sam_manyN_fname.append("_manyN.bam");

    std::ofstream out_agg(out_agg_fname);

    samFile *al = hts_open(args.get_string(Opt::INPUT).c_str(),"r");
    bam_hdr_t *al_hdr = sam_hdr_read(al);

    samFile *outSAM=sam_open(out_sam_fname.c_str(),"wb");
    bam_hdr_t *outSAM_header=bam_hdr_dup(al_hdr);
    int ret_hdr = sam_hdr_write(outSAM,outSAM_header);

    samFile *outSAM_manyN=sam_open(out_sam_manyN_fname.c_str(),"wb");
    bam_hdr_t *outSAM_manyN_header=bam_hdr_dup(al_hdr);
    ret_hdr = sam_hdr_write(outSAM_manyN,outSAM_manyN_header);

    bam1_t *curAl = bam_init1(); // initialize the alignment record

    while(sam_read1(al,al_hdr,curAl)>0) { // only perfom if unaligned flag is set to true
        if (curAl->core.flag & 4) { // if read is unmapped
            continue;
        }

        std::vector<int> coords;
        bool isSpliced = process_cigar(curAl,coords);
        if(isSpliced){
            bool found_wrong = false;
            for(int i=0;i<coords.size();i+=3){
                if(coords[i]<args.get_int(Opt::OVERHANG)||coords[i+2]<args.get_int(Opt::OVERHANG)){ // short overhang
                    found_wrong = true;
                    out_agg<<bam_get_qname(curAl)<<"\t"
                    <<al_hdr->target_name[curAl->core.tid]<<"\t"
                    <<curAl->core.pos<<"\t"
                    <<coords[i]<<"\t"
                    <<coords[i+1]<<"\t"
                    <<coords[i+2]<<std::endl;
                }
            }
            if(found_wrong){
                add_num_splice_tag(curAl,coords.size()/3);
                int ret_val = sam_write1(outSAM, outSAM_header, curAl);
                if(!ret_val){
                    std::cerr<<"something wrong when writing the read"<<std::endl;
                    exit(-1);
                }
            }
            if(coords.size()/3>=args.get_int(Opt::NUM_SPLICE)){ // too many splice sites
                add_num_splice_tag(curAl,coords.size()/3);
                int ret_val = sam_write1(outSAM_manyN, outSAM_manyN_header, curAl);
                if(!ret_val){
                    std::cerr<<"something wrong when writing the read"<<std::endl;
                    exit(-1);
                }
            }
        }
    }
    bam_destroy1(curAl);

    bam_hdr_destroy(al_hdr);
    sam_close(al);
    sam_close(outSAM);
    bam_hdr_destroy(outSAM_header);
    sam_close(outSAM_manyN);
    bam_hdr_destroy(outSAM_manyN_header);

    out_agg.close();

    return 0;
}