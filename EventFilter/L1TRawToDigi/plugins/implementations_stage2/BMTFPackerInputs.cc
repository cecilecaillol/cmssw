#include "BMTFPackerInputs.h"

#include <vector>
//#include <bitset>//debug

namespace l1t {
  namespace stage2 {

    const int BMTFPackerInputs::ownLinks_[] = {4, 5, 12, 13, 20, 21, 22, 23, 28, 29};

    Blocks BMTFPackerInputs::pack(const edm::Event& event, const PackerTokens* toks) {
      int board_id = (int)board();

      auto muonPhToken = static_cast<const BMTFTokens*>(toks)->getInputMuonTokenPh();
      auto muonThToken = static_cast<const BMTFTokens*>(toks)->getInputMuonTokenTh();

      Blocks blocks;

      edm::Handle<L1MuDTChambPhContainer> phInputs;
      event.getByToken(muonPhToken, phInputs);
      edm::Handle<L1MuDTChambThContainer> thInputs;
      event.getByToken(muonThToken, thInputs);

      uint32_t qualEta_32bit = 0;
      uint32_t posEta_0_32bit = 0;
      uint32_t posEta_n5_32bit = 0, posEta_n4_32bit = 0, posEta_n3_32bit = 0, posEta_n2_32bit = 0, posEta_n1_32bit = 0;
      uint32_t posEta_p5_32bit = 0, posEta_p4_32bit = 0, posEta_p3_32bit = 0, posEta_p2_32bit = 0, posEta_p1_32bit = 0;
      uint32_t posEta_n10_32bit = 0, posEta_n9_32bit = 0, posEta_n8_32bit = 0, posEta_n7_32bit = 0, posEta_n6_32bit = 0;
      uint32_t posEta_p10_32bit = 0, posEta_p9_32bit = 0, posEta_p8_32bit = 0, posEta_p7_32bit = 0, posEta_p6_32bit = 0;

      bool moreBXeta = false;
      for (int link = 0; link <= 35; link++) {
        //	  std::cout << "link : " << link << std::endl;

        if ((link >= 6 && link < 8) || (link >= 14 && link < 16) || (link >= 30 && link < 32))
          continue;

        //initializing null block_payloads and block_id
        std::vector<uint32_t> payload_0(6, 0);
        std::vector<uint32_t> payload_p1(6, 0);
        std::vector<uint32_t> payload_n1(6, 0);
        std::vector<uint32_t> payload_p2(6, 0);
        std::vector<uint32_t> payload_n2(6, 0);
	std::vector<uint32_t> payload_p3(6, 0);
        std::vector<uint32_t> payload_n3(6, 0);
	std::vector<uint32_t> payload_p4(6, 0);
        std::vector<uint32_t> payload_n4(6, 0);
	std::vector<uint32_t> payload_p5(6, 0);
        std::vector<uint32_t> payload_n5(6, 0);
	std::vector<uint32_t> payload_p6(6, 0);
        std::vector<uint32_t> payload_n6(6, 0);
        std::vector<uint32_t> payload_p7(6, 0);
        std::vector<uint32_t> payload_n7(6, 0);
        std::vector<uint32_t> payload_p8(6, 0);
        std::vector<uint32_t> payload_n8(6, 0);
        std::vector<uint32_t> payload_p9(6, 0);
        std::vector<uint32_t> payload_n9(6, 0);
        std::vector<uint32_t> payload_p10(6, 0);
        std::vector<uint32_t> payload_n10(6, 0);

        unsigned int block_id = 2 * link;

        std::vector<bool> bxPresent(21, false);
        bool moreBXphi = false;

        unsigned int BC = 0;

        //The first 4 phi words for the link's payload

        for (const auto& iphi : *(phInputs->getContainer())) {
          // Only allow -2 <= bxNum <= 2, as in Stage2 data
          if (std::abs(iphi.bxNum()) > 10) //FIXME Cecile
            continue;

          if (iphi.bxNum() != 0)
            moreBXphi = true;

          //	      std::cout << "scNum+1 = " << iphi->scNum()+1 << ",   board_id = " << board_id << std::endl;
          // 	      std::cout << "bx, station = " << iphi->bxNum() << ", " << iphi->stNum() << std::endl;

          //BC = iphi->BxCnt();//this thing here is not completely functional

          if (iphi.scNum() + 1 != board_id)
            continue;
          //	      std::cout << "correct board" << std::endl;

          if (link != ownLinks_[4 + 2 * (iphi.whNum()) + iphi.Ts2Tag()])
            continue;
          // 	      std::cout << "correct link" << std::endl;

          bxPresent[10 + iphi.bxNum()] = true;

          //1 create 32word, 2 insert 32word in correct Block Slot
          uint32_t word_32bit = wordPhMaker(iphi);  //1
	  if (bxPresent[0]) {
            payload_n10[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[1]) {
            payload_n9[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[2]) {
            payload_n8[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[3]) {
            payload_n7[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[4]) {
            payload_n6[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[5]) {
            payload_n5[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[6]) {
            payload_n4[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[7]) {
            payload_n3[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[8]) {
            payload_n2[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[9]) {
            payload_n1[iphi.stNum() - 1] = word_32bit;
	  } else if (bxPresent[10]) {
            payload_0[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[11]) {
            payload_p1[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[12]) {
            payload_p2[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[13]) {
            payload_p3[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[14]) {
            payload_p4[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[15]) {
            payload_p5[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[16]) {
            payload_p6[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[17]) {
            payload_p7[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[18]) {
            payload_p8[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[19]) {
            payload_p9[iphi.stNum() - 1] = word_32bit;
          } else if (bxPresent[20]) {
            payload_p10[iphi.stNum() - 1] = word_32bit;
          }

          bxPresent.assign(21, false);

        }  //phiCont_itr

        //============================================================================================

        //Create phiNull words (use the old format for null words because Packer gives fw_id=1)
        //uint32_t phiNull_32bit = 0 | (BC & 0x3) << 30 | (7 & 0x7) << 22; //for the new fw
        uint32_t phiNull_32bit = 0 | (BC & 0x3) << 30;  //for the old fw
        //phiNull = (BC)000001110000000000000000000000
        uint32_t etaNull_32bit = 0 | (BC & 0x3) << 30;
        //etaNull = (BC)000000000000000000000000000000

        //============================================================================================

        //The 5th & 6th words of the link's payload
        //	  std::cout << "link%2 = " << link%2 << std::endl;
        if (link % 2 == 0) {
          //these Eta vars have to be declared out of link itr scope in order to maintain its information for the next link
          //Using these as the basis for the pos and qual 32bit eta word
          //in case there are les than 3 hits, the entries will be zero.
          posEta_0_32bit = etaNull_32bit;
	  posEta_n5_32bit = etaNull_32bit;
	  posEta_n4_32bit = etaNull_32bit;
	  posEta_n3_32bit = etaNull_32bit;
          posEta_n2_32bit = etaNull_32bit;
          posEta_n1_32bit = etaNull_32bit;
	  posEta_p5_32bit = etaNull_32bit;
	  posEta_p4_32bit = etaNull_32bit;
	  posEta_p3_32bit = etaNull_32bit;
          posEta_p2_32bit = etaNull_32bit;
          posEta_p1_32bit = etaNull_32bit;

	  posEta_n10_32bit = etaNull_32bit;
          posEta_n9_32bit = etaNull_32bit;
          posEta_n8_32bit = etaNull_32bit;
          posEta_n7_32bit = etaNull_32bit;
          posEta_n6_32bit = etaNull_32bit;
          posEta_p10_32bit = etaNull_32bit;
          posEta_p9_32bit = etaNull_32bit;
          posEta_p8_32bit = etaNull_32bit;
          posEta_p7_32bit = etaNull_32bit;
          posEta_p6_32bit = etaNull_32bit;

          qualEta_32bit = etaNull_32bit;

          for (const auto& ithe : *(thInputs->getContainer())) {
            // Only allow -2 <= bxNum <= 2, as in Stage2 data
            if (std::abs(ithe.bxNum()) > 10)
              continue;

            if (ithe.bxNum() != 0)
              moreBXeta = true;

            //debug
            //		std::cout << "scNum+1 = " << ithe.scNum()+1 << ",   board_id = " << board_id << std::endl;
            //		std::cout << "bx, station = " << ithe.bxNum() << ", " << ithe.stNum() << std::endl;
            //		std::cout << "related link: " << ownLinks_[4+2*(ithe.whNum())] << std::endl;

            if (ithe.scNum() + 1 != board_id)
              continue;

            if (link != ownLinks_[4 + 2 * (ithe.whNum())])
              continue;

            bxPresent[10 + ithe.bxNum()] = true;

            //positions for the next link
            uint32_t posEta_7bit = wordThMaker(ithe, false);

            //qualities for this link
            uint32_t qualEta_7bit = wordThMaker(ithe, true);
            qualEta_32bit = qualEta_32bit | ((qualEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));

            //write the eta-pos and eta-qual information at the correct payload per BX
	    if (bxPresent[0]) {
              payload_n10[4] = qualEta_32bit;
              payload_n10[5] = etaNull_32bit | (2 & 0x2);
              posEta_n10_32bit = posEta_n10_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[1]) {
              payload_n9[4] = qualEta_32bit;
              payload_n9[5] = etaNull_32bit | (2 & 0x2);
              posEta_n9_32bit = posEta_n9_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[2]) {
              payload_n8[4] = qualEta_32bit;
              payload_n8[5] = etaNull_32bit | (2 & 0x2);
              posEta_n8_32bit = posEta_n8_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[3]) {
              payload_n7[4] = qualEta_32bit;
              payload_n7[5] = etaNull_32bit | (2 & 0x2);
              posEta_n7_32bit = posEta_n7_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[4]) {
              payload_n6[4] = qualEta_32bit;
              payload_n6[5] = etaNull_32bit | (2 & 0x2);
              posEta_n6_32bit = posEta_n6_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[5]) {
              payload_n5[4] = qualEta_32bit;
              payload_n5[5] = etaNull_32bit | (2 & 0x2);
              posEta_n5_32bit = posEta_n5_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[6]) {
              payload_n4[4] = qualEta_32bit;
              payload_n4[5] = etaNull_32bit | (2 & 0x2);
              posEta_n4_32bit = posEta_n4_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[7]) {
              payload_n3[4] = qualEta_32bit;
              payload_n3[5] = etaNull_32bit | (2 & 0x2);
              posEta_n3_32bit = posEta_n3_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[8]) {
              payload_n2[4] = qualEta_32bit;
              payload_n2[5] = etaNull_32bit | (2 & 0x2);
              posEta_n2_32bit = posEta_n2_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[9]) {
              payload_n1[4] = qualEta_32bit;
              payload_n1[5] = etaNull_32bit | (2 & 0x2);
              posEta_n1_32bit = posEta_n1_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[10]) {
              payload_0[4] = qualEta_32bit;
              payload_0[5] = etaNull_32bit | (2 & 0x2);
              posEta_0_32bit = posEta_0_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[11]) {
              payload_p1[4] = qualEta_32bit;
              payload_p1[5] = etaNull_32bit | (2 & 0x2);
              posEta_p1_32bit = posEta_p1_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[12]) {
              payload_p2[4] = qualEta_32bit;
              payload_p2[5] = etaNull_32bit | (2 & 0x2);
              posEta_p2_32bit = posEta_p2_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[13]) {
              payload_p3[4] = qualEta_32bit;
              payload_p3[5] = etaNull_32bit | (2 & 0x2);
              posEta_p3_32bit = posEta_p3_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[14]) {
              payload_p4[4] = qualEta_32bit;
              payload_p4[5] = etaNull_32bit | (2 & 0x2);
              posEta_p4_32bit = posEta_p4_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[15]) {
              payload_p5[4] = qualEta_32bit;
              payload_p5[5] = etaNull_32bit | (2 & 0x2);
              posEta_p5_32bit = posEta_p5_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[16]) {
              payload_p6[4] = qualEta_32bit;
              payload_p6[5] = etaNull_32bit | (2 & 0x2);
              posEta_p6_32bit = posEta_p6_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[17]) {
              payload_p7[4] = qualEta_32bit;
              payload_p7[5] = etaNull_32bit | (2 & 0x2);
              posEta_p7_32bit = posEta_p7_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[18]) {
              payload_p8[4] = qualEta_32bit;
              payload_p8[5] = etaNull_32bit | (2 & 0x2);
              posEta_p8_32bit = posEta_p8_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[19]) {
              payload_p9[4] = qualEta_32bit;
              payload_p9[5] = etaNull_32bit | (2 & 0x2);
              posEta_p9_32bit = posEta_p9_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            } else if (bxPresent[20]) {
              payload_p10[4] = qualEta_32bit;
              payload_p10[5] = etaNull_32bit | (2 & 0x2);
              posEta_p10_32bit = posEta_p10_32bit | ((posEta_7bit & 0x7F) << 7 * (ithe.stNum() - 1));
            }

            bxPresent.assign(21, false);

          }  //theCont_itr

        } else {  //now that we are in the next prime link #, write the buffered eta-qual

          if (moreBXeta) {
	    payload_n10[4] = posEta_n10_32bit;
            payload_n10[5] = etaNull_32bit;

	    payload_n9[4] = posEta_n9_32bit;
            payload_n9[5] = etaNull_32bit;

            payload_n8[4] = posEta_n8_32bit;
            payload_n8[5] = etaNull_32bit;

            payload_n7[4] = posEta_n7_32bit;
            payload_n7[5] = etaNull_32bit;

            payload_n6[4] = posEta_n6_32bit;
            payload_n6[5] = etaNull_32bit;

	    payload_n5[4] = posEta_n5_32bit;
            payload_n5[5] = etaNull_32bit;

            payload_n4[4] = posEta_n4_32bit;
            payload_n4[5] = etaNull_32bit;

            payload_n3[4] = posEta_n3_32bit;
            payload_n3[5] = etaNull_32bit;

            payload_n2[4] = posEta_n2_32bit;
            payload_n2[5] = etaNull_32bit;

            payload_n1[4] = posEta_n1_32bit;
            payload_n1[5] = etaNull_32bit;
          }

          payload_0[4] = posEta_0_32bit;
          payload_0[5] = etaNull_32bit;

          if (moreBXeta) {
            payload_p1[4] = posEta_p1_32bit;
            payload_p1[5] = etaNull_32bit;

            payload_p2[4] = posEta_p2_32bit;
            payload_p2[5] = etaNull_32bit;

	    payload_p3[4] = posEta_p3_32bit;
            payload_p3[5] = etaNull_32bit;

	    payload_p4[4] = posEta_p4_32bit;
            payload_p4[5] = etaNull_32bit;

	    payload_p5[4] = posEta_p5_32bit;
            payload_p5[5] = etaNull_32bit;

	    payload_p6[4] = posEta_p6_32bit;
            payload_p6[5] = etaNull_32bit;

            payload_p7[4] = posEta_p7_32bit;
            payload_p7[5] = etaNull_32bit;

            payload_p8[4] = posEta_p8_32bit;
            payload_p8[5] = etaNull_32bit;

            payload_p9[4] = posEta_p9_32bit;
            payload_p9[5] = etaNull_32bit;

            payload_p10[4] = posEta_p10_32bit;
            payload_p10[5] = etaNull_32bit;
          }
        }

        //std::cout << "moreBXphi" << moreBXphi << std::endl;
        //std::cout << "moreBXeta" << moreBXeta << std::endl;

        bool moreBX = moreBXphi || moreBXeta;
        //std::cout << "moreBX" << moreBX << std::endl;

        //case where phi words are notcreated
        for (int iSt = 0; iSt <= 3; iSt++) {

	  if (moreBX && payload_n10[iSt] == 0)
            payload_n10[iSt] = phiNull_32bit;

          if (moreBX && payload_n9[iSt] == 0)
            payload_n9[iSt] = phiNull_32bit;

          if (moreBX && payload_n8[iSt] == 0)
            payload_n8[iSt] = phiNull_32bit;

          if (moreBX && payload_n7[iSt] == 0)
            payload_n7[iSt] = phiNull_32bit;

          if (moreBX && payload_n6[iSt] == 0)
            payload_n6[iSt] = phiNull_32bit;

	  if (moreBX && payload_n5[iSt] == 0)
            payload_n5[iSt] = phiNull_32bit;

	  if (moreBX && payload_n4[iSt] == 0)
            payload_n4[iSt] = phiNull_32bit;

	  if (moreBX && payload_n3[iSt] == 0)
            payload_n3[iSt] = phiNull_32bit;

          if (moreBX && payload_n2[iSt] == 0)
            payload_n2[iSt] = phiNull_32bit;

          if (moreBX && payload_n1[iSt] == 0)
            payload_n1[iSt] = phiNull_32bit;

          if (payload_0[iSt] == 0)
            payload_0[iSt] = phiNull_32bit;

          if (moreBX && payload_p1[iSt] == 0)
            payload_p1[iSt] = phiNull_32bit;

          if (moreBX && payload_p2[iSt] == 0)
            payload_p2[iSt] = phiNull_32bit;

	  if (moreBX && payload_p3[iSt] == 0)
            payload_p3[iSt] = phiNull_32bit;

	  if (moreBX && payload_p4[iSt] == 0)
            payload_p4[iSt] = phiNull_32bit;

	  if (moreBX && payload_p5[iSt] == 0)
            payload_p5[iSt] = phiNull_32bit;

	  if (moreBX && payload_p6[iSt] == 0)
            payload_p6[iSt] = phiNull_32bit;

          if (moreBX && payload_p7[iSt] == 0)
            payload_p7[iSt] = phiNull_32bit;

          if (moreBX && payload_p8[iSt] == 0)
            payload_p8[iSt] = phiNull_32bit;

          if (moreBX && payload_p9[iSt] == 0)
            payload_p9[iSt] = phiNull_32bit;

          if (moreBX && payload_p10[iSt] == 0)
            payload_p10[iSt] = phiNull_32bit;
        }

        //case where eta words are not created
        for (int word = 4; word <= 5; word++) {
	  if (moreBX && payload_n10[word] == 0)
            payload_n10[word] = etaNull_32bit;

          if (moreBX && payload_n9[word] == 0)
            payload_n9[word] = etaNull_32bit;

          if (moreBX && payload_n8[word] == 0)
            payload_n8[word] = etaNull_32bit;

          if (moreBX && payload_n7[word] == 0)
            payload_n7[word] = etaNull_32bit;

          if (moreBX && payload_n6[word] == 0)
            payload_n6[word] = etaNull_32bit;

	  if (moreBX && payload_n5[word] == 0)
            payload_n5[word] = etaNull_32bit;

	  if (moreBX && payload_n4[word] == 0)
            payload_n4[word] = etaNull_32bit;

          if (moreBX && payload_n3[word] == 0)
            payload_n3[word] = etaNull_32bit;

          if (moreBX && payload_n2[word] == 0)
            payload_n2[word] = etaNull_32bit;

          if (moreBX && payload_n1[word] == 0)
            payload_n1[word] = etaNull_32bit;

          if (payload_0[word] == 0)
            payload_0[word] = etaNull_32bit;

          if (moreBX && payload_p1[word] == 0)
            payload_p1[word] = etaNull_32bit;

          if (moreBX && payload_p2[word] == 0)
            payload_p2[word] = etaNull_32bit;

	  if (moreBX && payload_p3[word] == 0)
            payload_p3[word] = etaNull_32bit;

	  if (moreBX && payload_p4[word] == 0)
            payload_p4[word] = etaNull_32bit;

	  if (moreBX && payload_p5[word] == 0)
            payload_p5[word] = etaNull_32bit;

	  if (moreBX && payload_p6[word] == 0)
            payload_p6[word] = etaNull_32bit;

          if (moreBX && payload_p7[word] == 0)
            payload_p7[word] = etaNull_32bit;

          if (moreBX && payload_p8[word] == 0)
            payload_p8[word] = etaNull_32bit;

          if (moreBX && payload_p9[word] == 0)
            payload_p9[word] = etaNull_32bit;

          if (moreBX && payload_p10[word] == 0)
            payload_p10[word] = etaNull_32bit;
        }

        //============================================================================================

        /*

	  //debug
	  std::cout << "payload created : " << std::endl;
	  if (moreBX){	    
	    for (auto &word : payload_n2)	// 
	      std::cout << std::bitset<32>(word).to_string() << std::endl;
	  
	    for (auto &word : payload_n1)	// 
	      std::cout << std::bitset<32>(word).to_string() << std::endl;
	  }
	  
	  for (auto &word : payload_0)	// 
	    std::cout << std::bitset<32>(word).to_string() << std::endl;
	  
	  if (moreBX){
	    for (auto &word : payload_p1)	// 
	      std::cout << std::bitset<32>(word).to_string() << std::endl;
	  
	    for (auto &word : payload_p2)	// 
	      std::cout << std::bitset<32>(word).to_string() << std::endl;
	  }
	  */
        //============================================================================================

        std::vector<uint32_t> payload;

        if (moreBX) {  //push -2,-1 bx payloads
	  for (int i = 0; i < 6; i++)
            payload.push_back(payload_n10[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n9[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n8[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n7[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n6[i]);
	  for (int i = 0; i < 6; i++)
            payload.push_back(payload_n5[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n4[i]);
	  for (int i = 0; i < 6; i++)
            payload.push_back(payload_n3[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n2[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_n1[i]);
        }

        for (int i = 0; i < 6; i++)  //zero bx payload
          payload.push_back(payload_0[i]);

        if (moreBX) {  //push +1,+2 bx payloads
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p1[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p2[i]);
	  for (int i = 0; i < 6; i++)
            payload.push_back(payload_p3[i]);
	  for (int i = 0; i < 6; i++)
            payload.push_back(payload_p4[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p5[i]);
	  for (int i = 0; i < 6; i++)
            payload.push_back(payload_p6[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p7[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p8[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p9[i]);
          for (int i = 0; i < 6; i++)
            payload.push_back(payload_p10[i]);
        }

        //in format Block(id,payload)
        blocks.push_back(Block(block_id, payload));

        if (link % 2 != 0) {
          moreBXeta = false;
        }

      }  //link_itr

      return blocks;
    }

    uint32_t BMTFPackerInputs::wordPhMaker(const L1MuDTChambPhDigi& phInput) {
      uint32_t temp(0);

      temp = (phInput.phi() & phiMask) << phiShift | (phInput.phiB() & phiBMask) << phiBShift |
             (phInput.code() & qualMask) << qualShift | (phInput.RpcBit() & rpcMask) << rpcShift | (0) << 29 |
             (phInput.BxCnt() & bxCntMask) << bxCntShift;

      //	std::cout<<"uint32_t temp is: "<<std::bitset<32>(temp).to_string()<<"---- phi="<<phInput.phi()<<", phiB()="<<phInput.phiB()<<", code(qual?)="<<phInput.code()<<", RPC="<<phInput.RpcBit()<<", BC="<<phInput.BxCnt()<<std::endl;
      return temp;
    }

    uint32_t BMTFPackerInputs::wordThMaker(const L1MuDTChambThDigi& thInput, const bool& qualFlag) {
      uint32_t temp(0);
      if (!qualFlag) {
        for (int i = 6; i >= 0; i--) {
          temp = temp << 1;
          temp = temp | (thInput.position(i) & 0x1);
        }
      } else {
        for (int i = 6; i >= 0; i--) {
          temp = temp << 1;
          temp = temp | (thInput.quality(i) & 0x1);
        }
      }

      return temp;
    }

  }  // namespace stage2
}  // namespace l1t
DEFINE_L1T_PACKER(l1t::stage2::BMTFPackerInputs);
