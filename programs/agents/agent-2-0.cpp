/********************************************************************
* Copyright (C) 2015 by Interstel Technologies, Inc.
*   and Hawaii Space Flight Laboratory.
*
* This file is part of the COSMOS/core that is the central
* module for COSMOS. For more information on COSMOS go to
* <http://cosmos-project.com>
*
* The COSMOS/core software is licenced under the
* GNU Lesser General Public License (LGPL) version 3 licence.
*
* You should have received a copy of the
* GNU Lesser General Public License
* If not, go to <http://www.gnu.org/licenses/>
*
* COSMOS/core is free software: you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public License
* as published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
*
* COSMOS/core is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
* Lesser General Public License for more details.
*
* Refer to the "licences" folder for further information on the
* condititons and terms to use this software.
********************************************************************/

/*! \file agent.cpp
* \brief Agent control program source
*/

//! \ingroup general
//! \defgroup agent_client Agent control program
//! This program allows communication with any of the Agents on the local network.
//! With it you can:
//! - list available Agents
//! - acquire the request list of specific Agents
//! - command specific Agents
//! - monitor Agent traffic

#include "support/configCosmos.h"
#include <stdlib.h>
#include "agent/agentclass.h"
#include "support/jsonlib.h"
#include "physics/physicslib.h"
#include "support/datalib.h"
#include "sys/stat.h"
#include "limits.h"
#include <iostream>

const int REQUEST_WAIT_TIME = 2;
const int SERVER_WAIT_TIME = 6;

string output;
string node_name = "";
string agent_name = "";


int main(int argc, char *argv[])
{
    int nbytes;
    beatstruc cbeat;
    Agent *agent;

    // dont' print debug messages
    //agent->debug_level = 0;

    agent = new Agent();
    if (agent->cinfo == nullptr)
    {
        fprintf(agent->get_debug_fd(), "%16.10f %s Failed to start Agent %s on Node %s Dated %s : %s\n",currentmjd(), mjd2iso8601(currentmjd()).c_str(), agent->getAgent().c_str(), agent->getNode().c_str(), utc2iso8601(data_ctime(argv[0])).c_str(), cosmos_error_string(NODE_ERROR_NODE).c_str());
        exit(NODE_ERROR_NODE);
    }


    // check command line arguments
    switch (argc)
    {
    case 1:
        {
            printf("Usage: agent [ list | dump [soh, beat, ###] | node_name agent_name \"request [ arguments ]\" ]\n");
            exit(1);
        }
        break;
    case 2:
    case 3:
    case 4:
    case 5:
        // agent dump request
        if (!strcmp(argv[1],"dump"))
        {
            double lmjd = 0., dmjd;
            string channel;
            Agent::AgentMessage message_type;
            Agent::messstruc message;
            int i;
            locstruc loc;

// JIMNOTE: this block will never be entered

            switch(argc)
            {
            case 5:
                agent_name = argv[4];
            case 4:
                node_name = argv[3];
            case 3:
                channel = argv[2];
                if (channel == "soh")
                {
                    message_type = Agent::AgentMessage::SOH;
                }
                else if (channel == "beat")
                {
                    message_type = Agent::AgentMessage::BEAT;
                }
                else if (channel == "request")
                {
                    message_type = Agent::AgentMessage::REQUEST;
                }
                else if (channel == "response")
                {
                    message_type = Agent::AgentMessage::RESPONSE;
                }
                else
                {
                    message_type = (Agent::AgentMessage)atoi(channel.c_str());
                }
                break;
            case 2:
                channel.clear();
                message_type = Agent::AgentMessage::ALL;
                break;
            }

            while (1)
            {
                int32_t iretn;
                if ((iretn=agent->readring(message, message_type, 1., Agent::Where::TAIL)) > 0)
                {
                    Agent::AgentMessage message_type_read = (Agent::AgentMessage)iretn;

                    // Skip if either not Agent::AgentMessage::ALL, or not desired AGENT_MESSAGE
                    if (!channel.empty() && message_type != message_type_read)
                    {
                        continue;
                    }

                    if (!node_name.empty() && node_name != message.meta.beat.node)
                    {
                        continue;
                    }

                    if (!agent_name.empty() && agent_name != message.meta.beat.proc)
                    {
                        continue;
                    }

                    switch (message_type_read)
                    {
                    case Agent::AgentMessage::SOH:
                        printf("[SOH]");
                        break;
                    case Agent::AgentMessage::BEAT:
                        printf("[BEAT]");
                        break;
                    case Agent::AgentMessage::REQUEST:
                        printf("[REQUEST]");
                        break;
                    case Agent::AgentMessage::RESPONSE:
                        printf("[RESPONSE]");
                        break;
                    default:
                        printf("[%d]",message_type_read);
                        break;
                    }

                    printf("%.15g:[%s:%s][%s:%u](%lu:%lu:%zu)\n",message.meta.beat.utc, message.meta.beat.node, message.meta.beat.proc, message.meta.beat.addr, message.meta.beat.port, message.jdata.size(), message.adata.size(), message.bdata.size());
                    printf("%s\n",message.jdata.c_str());
                    if (message_type_read < Agent::AgentMessage::BINARY)
                    {
                        if (!channel.empty())
                        {
                            printf("%s\n",message.adata.c_str());
                        }
                    }

                    if ((channel=="info") && message_type_read == Agent::AgentMessage::TRACK)
                    {
                        if (agent->cinfo->node.loc.utc > 0.)
                        {
                            if (lmjd > 0.)
                                dmjd = 86400.*(agent->cinfo->node.loc.utc-lmjd);
                            else
                                dmjd = 0.;
                            loc.pos.icrf.s = agent->cinfo->node.loc.pos.icrf.s;
                            loc.pos.utc = agent->cinfo->node.loc.utc;
                            pos_eci(&loc);
                            printf("%16.15g %6.4g %s %8.3f %8.3f %8.3f %5.1f %5.1f %5.1f\n",agent->cinfo->node.loc.utc,dmjd,agent->cinfo->node.name,DEGOF(loc.pos.geod.s.lon),DEGOF(loc.pos.geod.s.lat),loc.pos.geod.s.h,agent->cinfo->node.phys.powgen,agent->cinfo->node.phys.powuse,agent->cinfo->node.phys.battlev);
                            lmjd = agent->cinfo->node.loc.utc;
                        }
                    }

                    if ((channel=="imu") && message_type_read == Agent::AgentMessage::IMU) {
                        for (i=0; i<agent->cinfo->devspec.imu_cnt; i++) {
                            if (agent->cinfo->agent[0].beat.utc > 0.) {
                                if (lmjd > 0.)
                                    dmjd = 86400.*(agent->cinfo->agent[0].beat.utc-lmjd);
                                else
                                    dmjd = 0.;
                                printf("%.15g %.4g\n",loc.utc,dmjd);
                                lmjd = agent->cinfo->agent[0].beat.utc;
                            }
                        }
                    }
                }
                fflush(stdout);
            } //end infinite while loop
            break;
        }
        else if (!strcmp(argv[1],"list"))
        {
            size_t agent_count = 0;
            ElapsedTime et;
            agent->post(Agent::AgentMessage::REQUEST, "heartbeat");
            COSMOS_SLEEP(.5);
            do {
                if (agent->agent_list.size() > agent_count) {
                    for (size_t i=agent_count; i<agent->agent_list.size(); ++i) {
                        beatstruc cbeat = agent->agent_list[i];
						// NS1
                        //agent->send_request(cbeat,(char *)"getvalue {\"agent_pid\"}", output, REQUEST_WAIT_TIME);
						// NS2 -- works!
                        agent->send_request(cbeat, "get_value {\"agent[0].pid\"}", output, REQUEST_WAIT_TIME);
                        printf("[%lu] %.15g %s %s %s %hu %u\n",i,cbeat.utc,cbeat.node,cbeat.proc,cbeat.addr,cbeat.port,cbeat.bsz);
                        printf("\t%s\n",output.c_str());
                        fflush(stdout);
                    }
                    agent_count = agent->agent_list.size();
                }
                COSMOS_SLEEP(.1);
            } while (et.split() < SERVER_WAIT_TIME);
            exit(0);
            break;
        } else if (!strcmp(argv[1],"list_json")) {
            size_t agent_count = 0;
            ElapsedTime et;
            agent->post(Agent::AgentMessage::REQUEST, "heartbeat");
            COSMOS_SLEEP(.1);
            printf("{\"agent_list\":[");
            do {
                if (agent->agent_list.size() > agent_count) {
                    for (size_t i=agent_count; i<agent->agent_list.size(); ++i) {
                        beatstruc cbeat = agent->agent_list[i];
						// NS1
                        //agent->send_request(cbeat,(char *)"getvalue {\"agent_pid\"}", output, REQUEST_WAIT_TIME);
						// NS2 -- works!
                        agent->send_request(cbeat, "get_value {\"agent[0].pid\"}", output, REQUEST_WAIT_TIME);

                        if(i>0) printf(",");
                        printf("{\"agent_proc\": \"%s\", ", cbeat.proc);
                        printf("\"agent_utc\": %.15g, ", cbeat.utc);
                        printf("\"agent_node\": \"%s\", ", cbeat.node);
                        printf("\"agent_addr\": \"%s\", ", cbeat.addr);
                        printf("\"agent_port\": %hu, ", cbeat.port);
                        printf("\"agent_bsz\": %u, ", cbeat.bsz);
                        // HANDLE RESPONSE OUTPUT FORMAT
                        size_t status_pos;
                        if((status_pos= output.find("[OK]")  )!= string::npos){
                            if(output.at(0) == '{'){
                                if(status_pos - 1 >= 0 && output.at(status_pos - 1) == '}'){
                                    printf("\"output\": %s,", output.substr(0, status_pos).c_str());
                                } else {
                                    printf("\"output\": %s,", output.c_str());
                                }
                            } else {
                                printf("\"output\": \"%s\",", output.substr(status_pos ).c_str());
                            }
                            printf("\"status\": \"OK\"}");
                        } else if((status_pos = output.find("[NOK]") )!= string::npos){
                            printf("\"status\": \"NOK\"}");
                        } else {
                             printf("\"output\": %s }", output.c_str());
                        }
                        fflush(stdout);
                    }
                    fflush(stdout);
                    agent_count = agent->agent_list.size();
                }
                COSMOS_SLEEP(.1);
            } while (et.split() < SERVER_WAIT_TIME);
            printf("]}\n");
            exit(0);
            break;
	}

    default:
        if (!strcmp(argv[1],"dump")) {
            double lmjd = 0., dmjd;
            string channel;
            Agent::AgentMessage message_type;
            Agent::messstruc message;
            string header;
            int i;
            locstruc loc;

            if(argc == 3) {
                channel = argv[2];
                if (channel == "soh") {
                    message_type = Agent::AgentMessage::SOH;
                } else {
                    if (channel == "beat") {
                        message_type = Agent::AgentMessage::BEAT;
                    } else {
                        message_type = (Agent::AgentMessage)atoi(channel.c_str());
                    }
                }
            } else {
                channel.clear();
                message_type = Agent::AgentMessage::ALL;
            }

            while (1)
            {
                int32_t iretn;
                if ((iretn=agent->readring(message, Agent::AgentMessage::ALL, 1., Agent::Where::TAIL)) > 0) {
                    Agent::AgentMessage message_type_read = (Agent::AgentMessage)iretn;
                    // Skip if either not Agent::AgentMessage::ALL, or not desired AGENT_MESSAGE
                    if (!channel.empty() && message_type != message_type_read) { continue; }

                    header.resize(message.meta.jlength);
                    if (message_type_read < Agent::AgentMessage::BINARY)
                    {
						// NS1
                        memcpy(&header[0], message.adata.data(), message.meta.jlength);
                        json_clear_cosmosstruc(JSON_STRUCT_NODE, agent->cinfo);
                        json_clear_cosmosstruc(JSON_STRUCT_DEVICE, agent->cinfo);
                        json_parse(message.adata.c_str(), agent->cinfo);
						// NS2
                    } else {
                        memcpy(&header[0], message.bdata.data(), message.meta.jlength);
                    }

                    switch (message_type_read) {
                    	case Agent::AgentMessage::SOH:
                        	printf("[SOH]");
                        	break;
                    	case Agent::AgentMessage::BEAT:
                        	printf("[BEAT]");
                        	break;
                    	case Agent::AgentMessage::REQUEST:
                        	printf("[REQUEST]");
                        	break;
                    	case Agent::AgentMessage::RESPONSE:
                        	printf("[RESPONSE]");
                        	break;
                    	default:
                        	printf("[%d]",message_type_read);
                        	break;
 					}

                    printf("[%d] %.15g %s %s %s %hu %u\n",i,message.meta.beat.utc,message.meta.beat.node,message.meta.beat.proc,message.meta.beat.addr,message.meta.beat.port,message.meta.beat.bsz);

                    if (message_type_read < Agent::AgentMessage::BINARY && !channel.empty()) {
                        printf("%s\n",message.adata.c_str());
                    }

                    if ((channel=="info") && message_type_read == Agent::AgentMessage::TRACK) {
                        if (agent->cinfo->node.loc.utc > 0.) {
                            if (lmjd > 0.)
                                dmjd = 86400.*(agent->cinfo->node.loc.utc-lmjd);
                            else
                                dmjd = 0.;
                            loc.pos.icrf.s = agent->cinfo->node.loc.pos.icrf.s;
                            loc.pos.utc = agent->cinfo->node.loc.utc;
                            pos_eci(&loc);
                            printf("%16.15g %6.4g %s %8.3f %8.3f %8.3f %5.1f %5.1f %5.1f\n",agent->cinfo->node.loc.utc,dmjd,agent->cinfo->node.name,DEGOF(loc.pos.geod.s.lon),DEGOF(loc.pos.geod.s.lat),loc.pos.geod.s.h,agent->cinfo->node.phys.powgen,agent->cinfo->node.phys.powuse,agent->cinfo->node.phys.battlev);
                            lmjd = agent->cinfo->node.loc.utc;
                        }
                    }

                    if ((channel=="imu") && message_type_read == Agent::AgentMessage::IMU) {
                        for (i=0; i<agent->cinfo->devspec.imu_cnt; i++) {
                            if (agent->cinfo->agent[0].beat.utc > 0.) {
                                if (lmjd > 0.)
                                    dmjd = 86400.*(agent->cinfo->agent[0].beat.utc-lmjd);
                                else
                                    dmjd = 0.;
                                printf("%.15g %.4g\n",loc.utc,dmjd);
                                lmjd = agent->cinfo->agent[0].beat.utc;
                            }
                        }
                    }
                }
                fflush(stdout);
            } //end infinite while loop
        } else {

//        cbeat = agent->find_agent(argv[1], argv[2], SERVER_WAIT_TIME);
//        if (cbeat.exists)
        	if ((nbytes = agent->get_agent(argv[1], argv[2], SERVER_WAIT_TIME, cbeat)) > 0) {
            	if(argc == 3)
            	{
                	nbytes = agent->send_request(cbeat, "help", std::ref(output), REQUEST_WAIT_TIME);
                	printf("%s [%d]\n", output.c_str(), nbytes);
            	}
            	else
            	{
                	string request;
                	request = argv[3];
                	for (size_t i=0; i<(size_t)argc-4; ++i)
                	{
                    	request += " ";
                    	request += argv[i+4];
                	}
                	nbytes = agent->send_request(cbeat,request.c_str(), output, REQUEST_WAIT_TIME);
	//                printf("%s [%d]\n", output.c_str(), nbytes);
	//                printf("{\"request_output\": %s, \"bytes\": %d }\n", output.c_str(), nbytes);
                	// HANDLE RESPONSE OUTPUT FORMAT
                	printf("{");
                	size_t status_pos;
                	if((status_pos= output.find("[OK]")  )!= string::npos){
                    	if(output.at(0) == '{'){
                        	if(status_pos - 1 >= 0 && output.at(status_pos - 1) == '}'){
                            	printf("\"output\": %s,", output.substr(0, status_pos).c_str());
                        	} else {
                            	printf("\"output\": %s,", output.c_str());
                        	}
                    	} else {
                        	printf("\"output\": \"%s\",", output.substr(0,status_pos ).c_str());
                    	}
                    	printf("\"status\": \"OK\"}\n");
                	} else if((status_pos = output.find("[NOK]") )!= string::npos){
                    	printf("\"status\": \"NOK\"}\n");
                	} else {
                     	printf("\"output\": %s }\n", output.c_str());
                	}
            	}
        	} else {
            	if (!nbytes){
                	fprintf(stderr,"node-agent pair [%s:%s] not found\n",argv[1],argv[2]);
                	printf("{\"error\": \"node-agent pair [%s:%s] not found\" }\n",argv[1],argv[2]);
            	}
            	else
                	printf("Error: %d\n", nbytes);
        	}
    	}
    }
}