#ifndef BEACON_H
#define BEACON_H

#include <atomic>
#include <tuple>
#include "support/configCosmos.h"
#include "agent/agentclass.h"
#include "math/mathlib.h"
#include "support/enumlib.h"
//#include "support/transferclass.h"
#include "support/packetcomm.h"

// Alias of a function that updates a beacon struct
//  arg: string, beacon's namespace name
//  return: void
using update_beacon_func = void (*)(string);

/** 
 * Implement common functionalities for handling PacketComm packets
*/
namespace Cosmos {
    namespace Support {
        class Beacon
        {
        public:
            Beacon()
            {

            }

            enum class TypeId : uint8_t
                {
                None = 0,
                CPU1BeaconS = 10,
                CPU2BeaconS = 11,
                TsenBeaconS = 20,
                EPSCPUBeaconS = 30,
                EPSPVBeaconS = 31,
                EPSSWCHBeaconS = 32,
                EPSBATTBeaconS = 33,
                ADCSCPUBeaconS = 40,
                ADCSMTRBeaconS = 41,
                ADCSRWBeaconS = 42,
                ADCSmag_beacon = 43,
                ADCSGPSBeaconS = 44,
                ADCSSTTBeaconS = 45,
                ADCSSSENBeaconS = 46,
                ADCSSunBeaconS = 47,
                ADCSNadirBeaconS = 48,
                CPUBeaconL = 110,
                TsenBeaconL = 120,
                EPSBCREGBeaconL = 130,
                EPSSWCHBeaconL = 131,
                EPSBATTBeaconL = 133,
                ADCSMTRBeaconL = 140,
                ADCSRWBeaconL = 141,
                ADCSIMUBeaconL = 142,
                ADCSGPSBeaconL = 143,
                ADCSSTTBeaconL = 144,
                ADCSSSENBeaconL = 145,
                ADCSATTBeaconL = 146,
                };

            map<TypeId, string> TypeString = {
                {TypeId::CPU1BeaconS, "CPU1BeaconS"},
                {TypeId::CPU2BeaconS, "CPU2BeaconS"},
                {TypeId::TsenBeaconS, "TsenBeaconS"},
                {TypeId::EPSCPUBeaconS, "EPSCPUBeaconS"},
                {TypeId::EPSPVBeaconS, "EPSPVBeaconS"},
                {TypeId::EPSSWCHBeaconS, "EPSSWCHBeaconS"},
                {TypeId::EPSBATTBeaconS, "EPSBATTBeaconS"},
                {TypeId::ADCSCPUBeaconS, "ADCSCPUBeaconS"},
                {TypeId::ADCSMTRBeaconS, "ADCSMTRBeaconS"},
                {TypeId::ADCSRWBeaconS, "ADCSRWBeaconS"},
                {TypeId::ADCSmag_beacon, "ADCSmag_beacon"},
                {TypeId::ADCSGPSBeaconS, "ADCSGPSBeaconS"},
                {TypeId::ADCSSTTBeaconS, "ADCSSTTBeaconS"},
                {TypeId::ADCSSSENBeaconS, "ADCSSSENBeaconS"},
                {TypeId::ADCSSunBeaconS, "ADCSSunBeaconS"},
                {TypeId::ADCSNadirBeaconS, "ADCSNadirBeaconS"},
                {TypeId::CPUBeaconL, "CPUBeaconL"},
                {TypeId::TsenBeaconL, "TsenBeaconL"},
                {TypeId::EPSBCREGBeaconL, "EPSBCREGBeaconL"},
                {TypeId::EPSSWCHBeaconL, "EPSSWCHBeaconL"},
                {TypeId::EPSBATTBeaconL, "EPSBATTBeaconL"},
                {TypeId::ADCSMTRBeaconL, "ADCSMTRBeaconL"},
                {TypeId::ADCSRWBeaconL, "ADCSRWBeaconL"},
                {TypeId::ADCSGPSBeaconL, "ADCSGPSBeaconL"},
                {TypeId::ADCSIMUBeaconL, "ADCSIMUBeaconL"},
                {TypeId::ADCSSTTBeaconL, "ADCSSTTBeaconL"},
                {TypeId::ADCSSSENBeaconL, "ADCSSSENBeaconL"},
                {TypeId::ADCSATTBeaconL, "ADCSATTBeaconL"},
            };

            int32_t Init();
            int32_t EncodeBinary(TypeId type, cosmosstruc *cinfo, vector<uint8_t> &data);
            int32_t EncodeJson(TypeId type, cosmosstruc *cinfo, string& Contents);
            int32_t EncodeJson(TypeId type, cosmosstruc *cinfo, vector<uint8_t>& Contents);
            int32_t Decode(vector<uint8_t> &data, cosmosstruc *cinfo);
//            int32_t Encode(TypeId type);
//            int32_t JSONDecode(string& Contents);

            // Short Beacons
            struct __attribute__ ((packed)) cpu1_beacons
            {
                uint8_t type = 10;
                float met = 0.;
                float load = 0.;
                float memory = 0.;
                float disk = 0.;
            };

            struct __attribute__ ((packed)) cpu2_beacons
            {
                uint8_t type = 11;
                float met = 0.;
                uint32_t uptime = 0;
                uint32_t bootcount = 0;
                uint32_t initialdate = 0;
            };

            struct __attribute__ ((packed)) tsen_beacons
            {
                uint8_t type = 20;
                float met = 0.;
                float temp[3] = {0.};
            } ;

            struct __attribute__ ((packed)) epscpu_beacons
            {
                uint8_t type = 30;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) epsbcreg_beacons
            {
                uint8_t type = 31;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) epsswch_beacons
            {
                uint8_t type = 32;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) epsbatt_beacons
            {
                uint8_t type = 33;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) adcscpu_beacons
            {
                uint8_t type = 40;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) adcsmtr_beacons
            {
                uint8_t type = 41;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) adcsrw_beacons
            {
                uint8_t type = 42;
                float met = 0.;
                float omega1 = 0.;
                float omega2 = 0.;
                float omega3 = 0.;
            } ;

            struct __attribute__ ((packed)) adcsmag_beacon
            {
                uint8_t type = 43;
                float met = 0.;
                float magx = 0.;
                float magy = 0.;
                float magz = 0.;
            } ;
            //struct __attribute__ ((packed)) adcsimu_beacons
            //{
                //uint8_t type = 43;
                //float met = 0.;
                //float magx = 0.;
                //float magy = 0.;
                //float magz = 0.;
            //} ;

            struct __attribute__ ((packed)) adcsgps_beacons
            {
                uint8_t type = 44;
                float met = 0.;
                float geocx = 0.;
                float geocy = 0.;
                float geocz = 0.;
            } ;

            struct __attribute__ ((packed)) adcsstt_beacons
            {
                uint8_t type = 45;
                float met = 0.;
                float heading = 0.;
                float elevation = 0.;
                float bearing = 0.;
            } ;

            struct __attribute__ ((packed)) adcsssen_beacons
            {
                uint8_t type = 46;
                float met = 0.;
                float volt = 0.;
                float amp = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) adcssun_beacons
            {
                uint8_t type = 47;
                float met = 0.;
                float azimuth = 0.;
                float elevation = 0.;
                float temp = 0.;
            } ;

            struct __attribute__ ((packed)) adcsnadir_beacons
            {
                uint8_t type = 48;
                float met = 0.;
                float azimuth = 0.;
                float elevation = 0.;
                float temp = 0.;
            } ;

            // Long Beacons
            struct __attribute__ ((packed)) cpu_beacon
            {
                float uptime = 0.;
                uint32_t bootcount = 0;
                uint16_t mload = 0.;
                uint16_t mmemory = 0.;
                uint16_t mdisk = 0.;
                uint16_t ctemp = 0.;
            } ;

            static constexpr uint8_t cpu_count = 180 / sizeof(cpu_beacon);
            struct __attribute__ ((packed)) cpus_beacon
            {
                uint8_t type = 110;
                float met = 0.;
                uint32_t initialdate = 0;
                cpu_beacon cpu[cpu_count];
//                double last_updated;
//                double last_sent;
            };

            static constexpr uint8_t tsen_count = 180 / sizeof(uint16_t);
            struct __attribute__ ((packed)) tsen_beacon
            {
                uint8_t type = 120;
                float met = 0.;
                uint16_t ctemp[tsen_count] = {0};
            } ;

            struct __attribute__ ((packed)) epsbcreg_beacon
            {
                int16_t mvolt = 0.;
                int16_t mamp = 0.;
            } ;

            static constexpr uint8_t epsbcreg_count = 180 / sizeof(epsbcreg_beacon);
            struct __attribute__ ((packed)) epsbcregs_beacon
            {
                uint8_t type = 130;
                float met = 0.;
                epsbcreg_beacon bcreg[epsbcreg_count];
            };

            struct __attribute__ ((packed)) epsswch_beacon
            {
                int16_t mvolt = 0.;
                int16_t mamp = 0.;
            } ;

            static constexpr uint8_t epsswch_count = 180 / sizeof(epsswch_beacon);
            struct __attribute__ ((packed)) epsswchs_beacon
            {
                uint8_t type = 131;
                float met = 0.;
                epsswch_beacon swch[epsswch_count];
            };

//            struct __attribute__ ((packed)) epsswch2_beacon
//            {
//                uint8_t type = 132;
//                float met = 0.;
//                epsswch_beacon swch[16];
//            };

            struct __attribute__ ((packed)) epsbatt_beacon
            {
                int16_t mvolt = 0.;
                int16_t mamp = 0.;
                uint16_t cpercent = 0.;
                uint16_t ctemp = 0.;
            } ;

            static constexpr uint8_t epsbatt_count = 180 / sizeof(epsbatt_beacon);
            struct __attribute__ ((packed)) epsbatts_beacon
            {
                uint8_t type = 133;
                float met = 0.;
                epsbatt_beacon batt[epsbatt_count];
            };

            struct __attribute__ ((packed)) adcsmtr_beacon
            {
                float mom = 0.;
                float align[4] = {0.};
            };

            static constexpr uint8_t adcsmtr_count = 180 / sizeof(adcsmtr_beacon);
            struct __attribute__ ((packed)) adcsmtrs_beacon
            {
                uint8_t type = 140;
                float met = 0.;
                adcsmtr_beacon mtr[adcsmtr_count];
            };

            struct __attribute__ ((packed)) adcsrw_beacon
            {
                float omega = 0.;
                float alpha = 0.;
//                float moi[3] = {0.};
//                float align[4] = {0.};
            };

            static constexpr uint8_t adcsrw_count = 180 / sizeof(adcsrw_beacon);
            struct __attribute__ ((packed)) adcsrws_beacon
            {
                uint8_t type = 141;
                float met = 0.;
                adcsrw_beacon rw[adcsrw_count];
            };

            struct __attribute__ ((packed)) adcsimu_beacon
            {
                float theta[4] = {0.};
                float omega = 0.;
                float alpha = 0.;
                float accel = 0.;
                float bfield = 0.;
                float bdot = 0.;
                float align[4] = {0.};
            };

            static constexpr uint8_t adcsimu_count = 180 / sizeof(adcsimu_beacon);
            struct __attribute__ ((packed)) adcsimu1_beacon
            {
                uint8_t type = 142;
                float met = 0.;
                adcsimu_beacon imu[adcsimu_count];
            };

            struct __attribute__ ((packed)) adcsgps_beacon
            {
                double utc = 0.;
                double geoc[3] = {0.};
                float geocv[3] = {0.};
            };

            struct __attribute__ ((packed)) adcsgps1_beacon
            {
                uint8_t type = 143;
                float met = 0.;
                adcsgps_beacon gps[3];
            };

            struct __attribute__ ((packed)) adcsstt_beacon
            {
                float theta[4] = {0.};
                float omega[3] = {0.};
                float alpha[3] = {0.};
                float align[4] = {0.};
            };

            struct __attribute__ ((packed)) adcsstt1_beacon
            {
                uint8_t type = 144;
                float met = 0.;
                adcsstt_beacon stt[3];
            };

            struct __attribute__ ((packed)) adcsssen_beacon
            {
                float elevation = 0.;
                int8_t align = 0;
            };

            struct __attribute__ ((packed)) adcsssen1_beacon
            {
                uint8_t type = 145;
                float met = 0.;
                adcsssen_beacon ssen[10];
            };

            struct __attribute__ ((packed)) adcsatt_beacon
            {
                float vector[3] = {0.};
                float align[4] = {0.};
            };

            struct __attribute__ ((packed)) adcsatt1_beacon
            {
                uint8_t type = 146;
                float met = 0.;
                adcsatt_beacon sun;
                adcsatt_beacon earth;
                adcsssen_beacon coarse[10];
            };



            double get_interval();
            /// Register a new beacon
            int32_t add_beacon(const string& name, uint8_t type, size_t size);
            // Specify send pattern
            int32_t set_pattern(const vector<string>& pattern);
            /// Get next beacon in cycle of set beacon pattern as a SLIP-encoded PacketComm packet
            PacketComm get_next();
            /// Get name-specified beacon as SLIP-encoded PacketComm packet
            PacketComm get_from_name(const string& name);

            TypeId type;
            /// Data of interest
//            vector<uint8_t> data;

        private:
            // Reference to calling agent for accessing namespace
//            Agent *agent;
            // Function pointer to a function that will update the beacon(s)
            update_beacon_func update;
            // Time interval (in seconds) between beacon sends
//            std::atomic<double> interval;
            double interval;
            // Map beacon name to size of the beacon struct in bytes
            map<string, size_t> beacon_size;
            // Map beacon name to it type ID
            map<string, uint8_t> beacon_typeID;
            // Vector of predefined beacons
            vector<string> send_pattern;
            mutex send_pattern_mtx;
            int pattern_idx;
            int current_beacon;
//            cosmosstruc* cinfo;


            bool validate_short_beacon();
            bool validate_long_beacon();
        };
    }
}

#endif
