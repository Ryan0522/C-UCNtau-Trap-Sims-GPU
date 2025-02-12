/*  Written in 2017 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include "../inc/xorshift.h"
#include <stdint.h>
#include <string.h>

/* NOTE: as of 2017-10-08, this generator has a different multiplier (a
   fixed-point representation of the golden ratio), which eliminates
   linear dependencies from one of the lowest bits. The previous
   multiplier was 1181783497276652981 (M_8 in the paper). If you need to
   tell apart the two generators, you can refer to this generator as
   xorshift1024*φ and to the previous one as xorshift1024*M_8.

   This is a fast, high-quality generator. If 1024 bits of state are too
   much, try a xoroshiro128+ generator.

   Note that the two lowest bits of this generator are LFSRs of degree
   1024, and thus will fail binary rank tests. The other bits needs a much
   higher degree to be represented as LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

uint64_t s[16]; 
int p;

uint64_t next(void) {
	const uint64_t s0 = s[p];
	uint64_t s1 = s[p = (p + 1) & 15];
	s1 ^= s1 << 31; // a
	s[p] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30); // b,c
	return s[p] * 0x9e3779b97f4a7c13;
}


/* This is the jump function for the generator. It is equivalent
   to 2^512 calls to next(); it can be used to generate 2^512
   non-overlapping subsequences for parallel computations. */

void jump(void) {
	static const uint64_t JUMP[] = { 0x84242f96eca9c41d,
		0xa3c65b8776f96855, 0x5b34a39f070b5837, 0x4489affce4f31a1e,
		0x2ffeeb0a48316f40, 0xdc2d9891fe68c022, 0x3659132bb12fea70,
		0xaac17d8efa43cab8, 0xc4cb815590989b13, 0x5ee975283d71c93b,
		0x691548c86c1bd540, 0x7910c41d10a1e6a5, 0x0b5fc64563b3e2a8,
		0x047f7684e9fc949d, 0xb99181f2d8f685ca, 0x284600e3f30e38c3
	};

	uint64_t t[16] = { 0 };
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b)
				for(int j = 0; j < 16; j++)
					t[j] ^= s[(j + p) & 15];
			next();
		}

	for(int j = 0; j < 16; j++)
		s[(j + p) & 15] = t[j];
}

/* This seed comes from the NIST random beacon, supposedly a "truly random"
   source of bits.
BF1BD245CB65A1CA
34655F2977F91378
9F0BA21F0CE81DA2
E4AC0B168095E9A9
5EF65D0F70B24FFD
6EEB6C8A8AEA0190
C074669439A7F3EA
0A3E877CA70014B5
----------------
2B970F0E4ABC01E3
75B5BD001E1CF3E5
FAF82EE6AF09F7F3
59D381339A6960D1
B95CB26B6231CDD8
25F409BD914BB8AE
6C815344B3F27630
CBBDA012B11FF272
   Accessed from https://beacon.nist.gov/rest/record/1515692820 */
void initxorshift(int n) {
    switch(n) {
        case 0:
            s[0] = 0xBF1BD245CB65A1CA;
            s[1] = 0x34655F2977F91378;
            s[2] = 0x9F0BA21F0CE81DA2;
            s[3] = 0xE4AC0B168095E9A9;
            s[4] = 0x5EF65D0F70B24FFD;
            s[5] = 0x6EEB6C8A8AEA0190;
            s[6] = 0xC074669439A7F3EA;
            s[7] = 0x0A3E877CA70014B5;
            s[8] = 0x2B970F0E4ABC01E3;
            s[9] = 0x75B5BD001E1CF3E5;
            s[10] = 0xFAF82EE6AF09F7F3;
            s[11] = 0x59D381339A6960D1;
            s[12] = 0xB95CB26B6231CDD8;
            s[13] = 0x25F409BD914BB8AE;
            s[14] = 0x6C815344B3F27630;
            s[15] = 0xCBBDA012B11FF272;
            break;
        case 1: /*5A71830C90F2FD2FDA7C8AD74F0D0AD17D4DDBB4FBED6FCC04B104B7C12DA309E44BF2CB0CA4AD4CEE4E3D607B239D87026494F8B85559D36F666D2DD7DEB2D12AE86169720CCD257DF256183710FE09C2DFFD7169AA0C48B3E3C9322F3B7639424503B305C6E7F0ACD6AE4A5979CC79C4A368A85F95850AE955189CFF6B72CD*/
            s[0] = 0x5A71830C90F2FD2F;
            s[1] = 0xDA7C8AD74F0D0AD1;
            s[2] = 0x7D4DDBB4FBED6FCC;
            s[3] = 0x04B104B7C12DA309;
            s[4] = 0xE44BF2CB0CA4AD4C;
            s[5] = 0xEE4E3D607B239D87;
            s[6] = 0x026494F8B85559D3;
            s[7] = 0x6F666D2DD7DEB2D1;
            s[8] = 0x2AE86169720CCD25;
            s[9] = 0x7DF256183710FE09;
            s[10] = 0xC2DFFD7169AA0C48;
            s[11] = 0xB3E3C9322F3B7639;
            s[12] = 0x424503B305C6E7F0;
            s[13] = 0xACD6AE4A5979CC79;
            s[14] = 0xC4A368A85F95850A;
            s[15] = 0xE955189CFF6B72CD;
            break;
        case 2: /*BCE969A3656BAA70B3292C57012B99A5AFC0C9FF09FF040DB7555FACC50E8768E8CF093143D17B8A0C6FBF1AB61835A7CDDE75F11017A8202B5E2CF6A56F87C5F24CC454DEB834FA5F96489CB8D37452E7222BE44EB20113717B5CF3F741DE4DF000AB0928AB6602F2A245D04C400802AECA142F097FA287D8B64125930059FF*/
            s[0] = 0xBCE969A3656BAA70;
            s[1] = 0xB3292C57012B99A5;
            s[2] = 0xAFC0C9FF09FF040D;
            s[3] = 0xB7555FACC50E8768;
            s[4] = 0xE8CF093143D17B8A;
            s[5] = 0x0C6FBF1AB61835A7;
            s[6] = 0xCDDE75F11017A820;
            s[7] = 0x2B5E2CF6A56F87C5;
            s[8] = 0xF24CC454DEB834FA;
            s[9] = 0x5F96489CB8D37452;
            s[10] = 0xE7222BE44EB20113;
            s[11] = 0x717B5CF3F741DE4D;
            s[12] = 0xF000AB0928AB6602;
            s[13] = 0xF2A245D04C400802;
            s[14] = 0xAECA142F097FA287;
            s[15] = 0xD8B64125930059FF;
            break;
        case 3: /*3B59CFF9FD79088C768B1C33F4D36F534F36199D735DD5A47C2D8EB0AB8E69D6394E4DD40800C0A7CBC8B6CFAE51C5EF404BB3BA84C471CBCC2C2F799D789A8E04BD3497C1B28EA9D13B3F26F09656A6AE5840B06E2F44BB4B8D10EF7142624F6B01CDED7444C2D46A78331C406F5A777C0F7289E7EBC88AEEE2BFA4E588B52E*/
            s[0] = 0x3B59CFF9FD79088C;
            s[1] = 0x768B1C33F4D36F53;
            s[2] = 0x4F36199D735DD5A4;
            s[3] = 0x7C2D8EB0AB8E69D6;
            s[4] = 0x394E4DD40800C0A7;
            s[5] = 0xCBC8B6CFAE51C5EF;
            s[6] = 0x404BB3BA84C471CB;
            s[7] = 0xCC2C2F799D789A8E;
            s[8] = 0x04BD3497C1B28EA9;
            s[9] = 0xD13B3F26F09656A6;
            s[10] = 0xAE5840B06E2F44BB;
            s[11] = 0x4B8D10EF7142624F;
            s[12] = 0x6B01CDED7444C2D4;
            s[13] = 0x6A78331C406F5A77;
            s[14] = 0x7C0F7289E7EBC88A;
            s[15] = 0xEEE2BFA4E588B52E;
            break;
        case 4: /*B497976F2A6E7C0AA35B51132EE82741608713019BB1FE7D74AD44C3110917D8489F4BC22D93B47A9602463171D9521C7C36F565C9FDFC257B6ACA8254CC0CCE7511756EF8E97E9611FBBDC12985C2926F2D6F486C01A245D595BA74A581C226890E42DBAB360A5F51818C5AE7A3CC75BE21A01018B3849BE6C49F18158B5DB5*/
            s[0] = 0xB497976F2A6E7C0A;
            s[1] = 0xA35B51132EE82741;
            s[2] = 0x608713019BB1FE7D;
            s[3] = 0x74AD44C3110917D8;
            s[4] = 0x489F4BC22D93B47A;
            s[5] = 0x9602463171D9521C;
            s[6] = 0x7C36F565C9FDFC25;
            s[7] = 0x7B6ACA8254CC0CCE;
            s[8] = 0x7511756EF8E97E96;
            s[9] = 0x11FBBDC12985C292;
            s[10] = 0x6F2D6F486C01A245;
            s[11] = 0xD595BA74A581C226;
            s[12] = 0x890E42DBAB360A5F;
            s[13] = 0x51818C5AE7A3CC75;
            s[14] = 0xBE21A01018B3849B;
            s[15] = 0xE6C49F18158B5DB5;
            break;
        case 5: /*50A1C214F1CA4C20DC4454E4D9E74D9DB0FB9958C12BF6A8384CB1D9F307E89B8D44A3DDB327B8D3CA9F9BEFAC3C2277567655E63C9C46BA38C44D0F64895E3526BC9B0E864D681D137F3CC58CB8798F54CB6CEB83972A81B88AF4019DEFA29F05FB04D736A6EEBD82E71539D33D2A981E7ABF55E64CD441B95F1ACCA90FAE7D*/
            s[0] = 0x50A1C214F1CA4C20;
            s[1] = 0xDC4454E4D9E74D9D;
            s[2] = 0xB0FB9958C12BF6A8;
            s[3] = 0x384CB1D9F307E89B;
            s[4] = 0x8D44A3DDB327B8D3;
            s[5] = 0xCA9F9BEFAC3C2277;
            s[6] = 0x567655E63C9C46BA;
            s[7] = 0x38C44D0F64895E35;
            s[8] = 0x26BC9B0E864D681D;
            s[9] = 0x137F3CC58CB8798F;
            s[10] = 0x54CB6CEB83972A81;
            s[11] = 0xB88AF4019DEFA29F;
            s[12] = 0x05FB04D736A6EEBD;
            s[13] = 0x82E71539D33D2A98;
            s[14] = 0x1E7ABF55E64CD441;
            s[15] = 0xB95F1ACCA90FAE7D;
            break;
        case 6: /*26B32D1745AC3E41CDC3D3635EF51C8E34EF9B8BF1FEF343E33D63A334B08AAFEC3AA8D477FF7C38AAE0D4584D1188DFB5B9340DFD0F0039CE630E635ACA776145D6A0DCD73A8C19738EF32C98ACE5EBA6515346E495E111FDC72D9C31B2CB6BABE17EB5CF2F4B325952FE663C5011195128E0A94FC6F90DC937C67ED31F60DD*/
            s[0] = 0x26B32D1745AC3E41;
            s[1] = 0xCDC3D3635EF51C8E;
            s[2] = 0x34EF9B8BF1FEF343;
            s[3] = 0xE33D63A334B08AAF;
            s[4] = 0xEC3AA8D477FF7C38;
            s[5] = 0xAAE0D4584D1188DF;
            s[6] = 0xB5B9340DFD0F0039;
            s[7] = 0xCE630E635ACA7761;
            s[8] = 0x45D6A0DCD73A8C19;
            s[9] = 0x738EF32C98ACE5EB;
            s[10] = 0xA6515346E495E111;
            s[11] = 0xFDC72D9C31B2CB6B;
            s[12] = 0xABE17EB5CF2F4B32;
            s[13] = 0x5952FE663C501119;
            s[14] = 0x5128E0A94FC6F90D;
            s[15] = 0xC937C67ED31F60DD;
            break;
        case 7: /*3AA8C7E10A60049DEA6779CFE2F3A51C112AC4D25FC55302243BF6972C9851CF9D89B4CFD993CE504B9DDDBE1F9AA400B5B00730483D187EEAC71790B0451B72DA72FEF61E18DE4DCA8C9F459B1362859B7725B4982E5AE101EC3A5956A33D0CC12FC23E25705771CF9FBC8EBB4FC84F815D35510CAE52F8C0D9A8754AC76FB0*/
            s[0] = 0x3AA8C7E10A60049D;
            s[1] = 0xEA6779CFE2F3A51C;
            s[2] = 0x112AC4D25FC55302;
            s[3] = 0x243BF6972C9851CF;
            s[4] = 0x9D89B4CFD993CE50;
            s[5] = 0x4B9DDDBE1F9AA400;
            s[6] = 0xB5B00730483D187E;
            s[7] = 0xEAC71790B0451B72;
            s[8] = 0xDA72FEF61E18DE4D;
            s[9] = 0xCA8C9F459B136285;
            s[10] = 0x9B7725B4982E5AE1;
            s[11] = 0x01EC3A5956A33D0C;
            s[12] = 0xC12FC23E25705771;
            s[13] = 0xCF9FBC8EBB4FC84F;
            s[14] = 0x815D35510CAE52F8;
            s[15] = 0xC0D9A8754AC76FB0;
            break;
        case 8: /*89D078FBC8CB185F59ECC1B7567F1C1E23879808B35E56134518A200B19B8AF7638BD7A0CF64FC0DC82C63627DB3E814D7EA5168C15EDB858E2E286CE7D86FE9CDB7D815F9A1B2D232BC584A7031EC48BDD50C45CFADA3369EFA3BCBA705DAF0020BD76639F8C1CA193FD8ED1E376BA922D067575A05405594AEDCF42B80C48E*/
            s[0] = 0x89D078FBC8CB185F;
            s[1] = 0x59ECC1B7567F1C1E;
            s[2] = 0x23879808B35E5613;
            s[3] = 0x4518A200B19B8AF7;
            s[4] = 0x638BD7A0CF64FC0D;
            s[5] = 0xC82C63627DB3E814;
            s[6] = 0xD7EA5168C15EDB85;
            s[7] = 0x8E2E286CE7D86FE9;
            s[8] = 0xCDB7D815F9A1B2D2;
            s[9] = 0x32BC584A7031EC48;
            s[10] = 0xBDD50C45CFADA336;
            s[11] = 0x9EFA3BCBA705DAF0;
            s[12] = 0x020BD76639F8C1CA;
            s[13] = 0x193FD8ED1E376BA9;
            s[14] = 0x22D067575A054055;
            s[15] = 0x94AEDCF42B80C48E;
            break;
        case 9: /*7BA61DD6DD4046A6ACE0550BC45EAE470CD1B53E50EC5DA02D916CB0B2BEE21AAA8720AE6E386B81201E83D2EA524196AE67D525151DD6D56C577850C218E26E246EFD758580B2E8506F9E4A04772E2C000BF1A7CB607D3FE0177A799967B48110FCF9AAB662B2A1BAC3221C0FD69512AE66AC294326C2E0CFB1EEB09A0D2E61*/
            s[0] = 0x7BA61DD6DD4046A6;
            s[1] = 0xACE0550BC45EAE47;
            s[2] = 0x0CD1B53E50EC5DA0;
            s[3] = 0x2D916CB0B2BEE21A;
            s[4] = 0xAA8720AE6E386B81;
            s[5] = 0x201E83D2EA524196;
            s[6] = 0xAE67D525151DD6D5;
            s[7] = 0x6C577850C218E26E;
            s[8] = 0x246EFD758580B2E8;
            s[9] = 0x506F9E4A04772E2C;
            s[10] = 0x000BF1A7CB607D3F;
            s[11] = 0xE0177A799967B481;
            s[12] = 0x10FCF9AAB662B2A1;
            s[13] = 0xBAC3221C0FD69512;
            s[14] = 0xAE66AC294326C2E0;
            s[15] = 0xCFB1EEB09A0D2E61;
            break;
        case 10: /*541EF5B4AAB2DE7394F92F1A7DF6AA5EFA2E48882B39B772FAFF92DBAE04D4686C4A0B9A02BED91C08787F8DEE7C01E1B7E438887D4EA7187D5462EC9816F04AB1446D6D7E17B4F0FBD0A454A1864DFA09AE65F449A5AF87D5138298B6499F727C5AD07DAD419202E51369F821A2092057E8AD730870F087197527D19A03B9F0*/
            s[0] = 0x541EF5B4AAB2DE73;
            s[1] = 0x94F92F1A7DF6AA5E;
            s[2] = 0xFA2E48882B39B772;
            s[3] = 0xFAFF92DBAE04D468;
            s[4] = 0x6C4A0B9A02BED91C;
            s[5] = 0x08787F8DEE7C01E1;
            s[6] = 0xB7E438887D4EA718;
            s[7] = 0x7D5462EC9816F04A;
            s[8] = 0xB1446D6D7E17B4F0;
            s[9] = 0xFBD0A454A1864DFA;
            s[10] = 0x09AE65F449A5AF87;
            s[11] = 0xD5138298B6499F72;
            s[12] = 0x7C5AD07DAD419202;
            s[13] = 0xE51369F821A20920;
            s[14] = 0x57E8AD730870F087;
            s[15] = 0x197527D19A03B9F0;
            break;
        case 11: /*9583315E355DD2D2C0CFBADAC22EA7E7739B6A9D8AD6FCFB388E582E208FA168EFEF69B652D7394558848791B2D63707BB327D640C2C1FF32E93C646EC83C8FABB3AA9168F052AEA72AED158AC4481E3ABD420F72C13B8CC20A450E87BD2F75D544980185E9ED6AD274A946F05B8B1A9145BABA7BE2BC05C9FC41A5F30A8EDDD*/
            s[0] = 0x9583315E355DD2D2;
            s[1] = 0xC0CFBADAC22EA7E7;
            s[2] = 0x739B6A9D8AD6FCFB;
            s[3] = 0x388E582E208FA168;
            s[4] = 0xEFEF69B652D73945;
            s[5] = 0x58848791B2D63707;
            s[6] = 0xBB327D640C2C1FF3;
            s[7] = 0x2E93C646EC83C8FA;
            s[8] = 0xBB3AA9168F052AEA;
            s[9] = 0x72AED158AC4481E3;
            s[10] = 0xABD420F72C13B8CC;
            s[11] = 0x20A450E87BD2F75D;
            s[12] = 0x544980185E9ED6AD;
            s[13] = 0x274A946F05B8B1A9;
            s[14] = 0x145BABA7BE2BC05C;
            s[15] = 0x9FC41A5F30A8EDDD;
            break;
        case 12: /*CAB048D3CE5700B161CAD2954DCA6D1A799BB2A9BE9165C2E1014EE07DDA33AAE245BA0D5BA2B30E27A3911FAE14787FC619B52620E066651B9C833A4D42AE8EDDB0F8EFD96E9D43AC468AB3FECC6597B48F9EDBEFC3D57A6AF7AB3F294F39904224CE677D273115196BBACC2CBE8A0562BB9D9CD448215FF10151025769EF08*/
            s[0] = 0xCAB048D3CE5700B1;
            s[1] = 0x61CAD2954DCA6D1A;
            s[2] = 0x799BB2A9BE9165C2;
            s[3] = 0xE1014EE07DDA33AA;
            s[4] = 0xE245BA0D5BA2B30E;
            s[5] = 0x27A3911FAE14787F;
            s[6] = 0xC619B52620E06665;
            s[7] = 0x1B9C833A4D42AE8E;
            s[8] = 0xDDB0F8EFD96E9D43;
            s[9] = 0xAC468AB3FECC6597;
            s[10] = 0xB48F9EDBEFC3D57A;
            s[11] = 0x6AF7AB3F294F3990;
            s[12] = 0x4224CE677D273115;
            s[13] = 0x196BBACC2CBE8A05;
            s[14] = 0x62BB9D9CD448215F;
            s[15] = 0xF10151025769EF08;
            break;
        default:
            s[0] = 0;
            s[1] = 0;
            s[2] = 0;
            s[3] = 0;
            s[4] = 0;
            s[5] = 0;
            s[6] = 0;
            s[7] = 0;
            s[8] = 0;
            s[9] = 0;
            s[10] = 0;
            s[11] = 0;
            s[12] = 0;
            s[13] = 0;
            s[14] = 0;
            s[15] = 0;
            break;
    }
}

/* This function returns numbers on [0.1), via http://xoroshiro.di.unimi.it/*/
double nextU01() {
    uint64_t u = next();
    return (u >> 11) * (1. / (UINT64_C(1) << 53));
}