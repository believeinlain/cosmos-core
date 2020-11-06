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

#include "device/general/usrp_lib.h"
#include "support/jsondef.h"

int32_t usrp_connect(string device, uint8_t address, usrp_handle &handle)
{
    int32_t iretn;

    handle.serial = new Serial(device.c_str(), USRP_BAUD, USRP_BITS, USRP_PARITY, USRP_STOPBITS);

    if (!handle.serial->get_open())
    {
        return (SERIAL_ERROR_OPEN);
    }

    handle.address = address;
    iretn = usrp_check_address(handle);

    return iretn;
}

int32_t usrp_disconnect(usrp_handle &handle)
{
    if (!handle.serial->get_open()) return (SERIAL_ERROR_OPEN);

    handle.serial->close_device();
    return 0;
}

int32_t usrp_write_header(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = handle.serial->put_char(0xfe);
    if (iretn < 0)
    {
        return iretn;
    }
    iretn = handle.serial->put_char(0xfe);
    if (iretn < 0)
    {
        return iretn;
    }
    iretn = handle.serial->put_char(handle.address);
    if (iretn < 0)
    {
        return iretn;
    }
    iretn = handle.serial->put_char(0xe0);
    if (iretn < 0)
    {
        return iretn;
    }

    return 0;
}

int32_t usrp_write(usrp_handle &handle, uint8_t command)
{
    int32_t iretn = 0;

    vector <uint8_t> data;
    iretn = usrp_write(handle, command, data);
    return iretn;
}

int32_t usrp_write(usrp_handle &handle, uint8_t command, uint8_t subcommand)
{
    int32_t iretn = 0;

    vector <uint8_t> data;
    iretn = usrp_write(handle, command, subcommand, data);
    return iretn;
}

int32_t usrp_write(usrp_handle &handle, uint8_t command, vector <uint8_t> message)
{
    int32_t iretn = 0;

    handle.mut.lock();

    iretn = usrp_write_header(handle);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }

    iretn = handle.serial->put_char(command);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }

    if (message.size())
    {
        iretn = handle.serial->put_data(static_cast <uint8_t *>(message.data()), message.size());
        if (iretn < 0)
        {
            handle.mut.unlock();
            return iretn;
        }
    }

    iretn = handle.serial->put_char(0xfd);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }

    uint8_t buffer[100];
    iretn = handle.serial->get_data( buffer, 100);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }
    size_t base = message.size() + 6;

    if (static_cast <size_t>(base) < base)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_OUTPUT;
    }

    if (static_cast <size_t>(base) == base)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_MISMATCH;
    }

    if (buffer[base] != 0xfe || buffer[base+1] != 0xfe || buffer[base+2] != 0xe0 || buffer[base+3] != handle.address || buffer[static_cast <size_t>(base)-1] != 0xfd)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_OUTPUT;
    }

    if (buffer[base+4] == 0xfa)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_MISMATCH;
    }

    if (buffer[base+4] == 0xfb)
    {
        handle.response.resize(0);
        handle.mut.unlock();
        return 0;
    }
    else
    {
        base += 5;
        handle.response.resize(static_cast <size_t>(base)-(base+1));
        memcpy(static_cast<void *>(handle.response.data()), &buffer[base], static_cast <size_t>(base)-(base+1));
        handle.mut.unlock();
        return static_cast <int32_t>(base-(base+1));
    }

}

int32_t usrp_write(usrp_handle &handle, uint8_t command, uint8_t subcommand, vector <uint8_t> message)
{
    int32_t iretn = 0;

    handle.mut.lock();

    iretn = usrp_write_header(handle);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }

    iretn = handle.serial->put_char(command);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }
    iretn = handle.serial->put_char(subcommand);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }
    if (message.size())
    {
        iretn = handle.serial->put_data(message);
        if (iretn < 0)
        {
            handle.mut.unlock();
            return iretn;
        }
    }
    iretn = handle.serial->put_char(0xfd);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }

    vector <uint8_t> buffer;
    iretn = handle.serial->get_data(buffer, 100);
    if (iretn < 0)
    {
        handle.mut.unlock();
        return iretn;
    }

    size_t base = message.size() + 7;
    if (static_cast <size_t>(iretn) < base)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_OUTPUT;
    }

    if (static_cast <size_t>(iretn) == base)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_MISMATCH;
    }

    if (buffer[base] != 0xfe || buffer[base+1] != 0xfe || buffer[base+2] != 0xe0 || buffer[base+3] != handle.address || buffer[static_cast <size_t>(base)-1] != 0xfd)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_OUTPUT;
    }

    if (buffer[base+4] == 0xfa)
    {
        handle.mut.unlock();
        return GENERAL_ERROR_MISMATCH;
    }

    if (buffer[base+4] == 0xfb)
    {
        handle.response.resize(0);
        handle.mut.unlock();
        return 0;
    }
    else
    {
        base += 6;
        handle.response.resize(static_cast <size_t>(base)-(base+1));
        memcpy(static_cast <void *>(handle.response.data()), &buffer[base], static_cast <size_t>(base)-(base+1));
        handle.mut.unlock();
        return static_cast <int32_t>((base)-(base+1));
    }
}

uint8_t usrp_byte(vector <uint8_t> response)
{
    uint8_t result = 0.;
    for (size_t i=0; i<2; ++i)
    {
        result *= 100.;
        result += 10. * (response[i] >> 4) + (response[i] % 16);
    }

    return result;
}

uint8_t usrp_freq2band(double frequency)
{
    uint8_t freqband;

    if (frequency < 1.8e6)
    {
        freqband = 14;
    }
    else if (frequency < 2.0e6)
    {
        freqband = 1;
    }
    else if (frequency >= 3.4e6 && frequency < 4.1e6)
    {
        freqband = 2;
    }
    else if (frequency >= 6.9e6 && frequency < 7.5e6)
    {
        freqband = 3;
    }
    else if (frequency >= 9.9e6 && frequency < 10.5e6)
    {
        freqband = 4;
    }
    else if (frequency >= 13.9e6 && frequency < 14.5e6)
    {
        freqband = 5;
    }
    else if (frequency >= 17.9e6 && frequency < 18.5e6)
    {
        freqband = 6;
    }
    else if (frequency >= 20.9e6 && frequency < 21.5e6)
    {
        freqband = 7;
    }
    else if (frequency >= 24.4e6 && frequency < 25.1e6)
    {
        freqband = 8;
    }
    else if (frequency >= 28.0e6 && frequency < 30.0e6)
    {
        freqband = 9;
    }
    else if (frequency >= 50.0e6 && frequency <= 54.0e6)
    {
        freqband = 10;
    }
    else if (frequency >= 108.0e6 && frequency <= 174.0e6)
    {
        freqband = 11;
    }
    else if (frequency >= 420.0e6 && frequency <= 480.0e6)
    {
        freqband = 12;
    }
    else if (frequency >= 1240.0e6 && frequency <1320.0e6)
    {
        freqband = 13;
    }
    else
    {
        freqband = 14;
    }
    return freqband;
}

int32_t usrp_set_bandpass(usrp_handle &handle, double bandpass)
{
    int32_t iretn = 0;

    iretn = usrp_get_mode(handle);

    if (iretn < 0)
    {
        return iretn;
    }

    uint8_t filtband = 1;

    switch (handle.mode)
    {
    case USRP_MODE_AM:
        if (bandpass < 200 || bandpass > 10000)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        if (bandpass < 3333.)
        {
            filtband = 3;
        }
        else if (bandpass < 6666.)
        {
            filtband = 2;
        }
        else
        {
            filtband = 1;
        }
        break;
    case USRP_MODE_CW:
    case USRP_MODE_CWR:
        if (bandpass < 50 || bandpass > 3600)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        if (bandpass < 300.)
        {
            filtband = 3;
        }
        else if (bandpass < 600.)
        {
            filtband = 2;
        }
        else
        {
            filtband = 1;
        }
        break;
    case USRP_MODE_DV:
    case USRP_MODE_FM:
        if (bandpass != 7000. && bandpass != 10000. && bandpass != 15000.)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        if (bandpass == 7000.)
        {
            filtband = 3;
        }
        else if (bandpass == 10000.)
        {
            filtband = 2;
        }
        else
        {
            filtband = 1;
        }
        break;
    case USRP_MODE_LSB:
    case USRP_MODE_USB:
        if (bandpass < 50 || bandpass > 3600)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        if (!handle.datamode)
        {
            if (bandpass < 2100.)
            {
                filtband = 3;
            }
            else if (bandpass < 2700.)
            {
                filtband = 2;
            }
            else
            {
                filtband = 1;
            }
        }
        else
        {
            if (bandpass < 300.)
            {
                filtband = 3;
            }
            else if (bandpass < 600.)
            {
                filtband = 2;
            }
            else
            {
                filtband = 1;
            }
        }
        break;
    case USRP_MODE_RTTY:
    case USRP_MODE_RTTYR:
        if (bandpass < 50 || bandpass > 2700)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        if (bandpass < 300.)
        {
            filtband = 2;
        }
        else if (bandpass < 600.)
        {
            filtband = 2;
        }
        else
        {
            filtband = 1;
        }
        break;
    default:
        return GENERAL_ERROR_OUTOFRANGE;
    }

    if (bandpass >= 1e10 || bandpass < 0)
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }

    iretn = usrp_set_mode(handle, handle.opmode, filtband);
    if (iretn < 0)
    {
        return iretn;
    }

    handle.bandpass = static_cast <float>(bandpass);
    handle.filtband = filtband;

    return 0;
}

int32_t usrp_get_bandpass(usrp_handle &handle)
{
    int32_t iretn;

    iretn = usrp_get_mode(handle);
    if (iretn < 0)
    {
        return iretn;
    }

    switch (handle.mode)
    {
    case USRP_MODE_AM:
        switch(handle.filtband)
        {
        case 3:
            handle.bandpass = 3000.;
            break;
        case 2:
            handle.bandpass = 6000.;
            break;
        case 1:
            handle.bandpass = 9000.;
            break;
        }
        break;
    case USRP_MODE_CW:
    case USRP_MODE_CWR:
        switch(handle.filtband)
        {
        case 3:
            handle.bandpass = 250.;
            break;
        case 2:
            handle.bandpass = 500.;
            break;
        case 1:
            handle.bandpass = 1200.;
            break;
        }
        break;
    case USRP_MODE_DV:
    case USRP_MODE_FM:
        switch(handle.filtband)
        {
        case 3:
            handle.bandpass = 7000.;
            break;
        case 2:
            handle.bandpass = 10000.;
            break;
        case 1:
            handle.bandpass = 15000.;
            break;
        }
        break;
    case USRP_MODE_LSB:
    case USRP_MODE_USB:
        if (!handle.datamode)
        {
            switch(handle.filtband)
            {
            case 3:
                handle.bandpass = 1800.;
                break;
            case 2:
                handle.bandpass = 2400.;
                break;
            case 1:
                handle.bandpass = 3000.;
                break;
            }
        }
        else
        {
            switch(handle.filtband)
            {
            case 3:
                handle.bandpass = 350.;
                break;
            case 2:
                handle.bandpass = 500.;
                break;
            case 1:
                handle.bandpass = 1200.;
                break;
            }
        }
        break;
    case USRP_MODE_RTTY:
    case USRP_MODE_RTTYR:
        switch(handle.filtband)
        {
        case 3:
            handle.bandpass = 250.;
            break;
        case 2:
            handle.bandpass = 500.;
            break;
        case 1:
            handle.bandpass = 2400.;
            break;
        }
        break;
    }

    return 0;
}

int32_t usrp_set_rfgain(usrp_handle &handle, uint8_t rfgain)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    vector <uint8_t> data { 0x0,0x0 };

    for (size_t i=0; i<2; ++i)
    {
        data[1-i] = 0;
        for (size_t j=0; j<2; ++j)
        {
            uint8_t digit = rfgain % 10;
            switch (j)
            {
            case 0:
                data[1-i] += digit;
                break;
            case 1:
                data[1-i] += digit << 4;
                break;
            }
            rfgain /= 10;
        }
    }

    iretn = usrp_write(handle, 0x14, 0x2, data);
    return iretn;
}

int32_t usrp_set_squelch(usrp_handle &handle, uint8_t squelch)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    vector <uint8_t> data { 0x0,0x0 };

    for (size_t i=0; i<2; ++i)
    {
        data[1-i] = 0;
        for (size_t j=0; j<2; ++j)
        {
            uint8_t digit = squelch % 10;
            switch (j)
            {
            case 0:
                data[1-i] += digit;
                break;
            case 1:
                data[1-i] += digit << 4;
                break;
            }
            squelch /= 10;
        }
    }

    iretn = usrp_write(handle, 0x14, 0x3, data);
    return iretn;
}

int32_t usrp_set_rfpower(usrp_handle &handle, float power)
{
    uint8_t rfpower = 0;
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    if (handle.freqband < 11)
    {
        if (power < 2.f || power > (handle.mode==USRP_MODE_AM?30.f:100.f))
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        rfpower = static_cast <uint8_t> (255 * (power - 2.f) / (handle.mode==USRP_MODE_AM?28.f:98.f));
    }
    else if (handle.freqband < 12)
    {
        if (power < 2.f || power > (100.f))
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        rfpower = static_cast <uint8_t> (255 * (power - 2.f) / (98.f));
    }
    else if (handle.freqband < 13)
    {
        if (power < 2.f || power > 75.f)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        rfpower = static_cast <uint8_t> (255 * (power - 2.f) / 73.f);
    }
    else if (handle.freqband < 14)
    {
        if (power < 2.f || power > 10.f)
        {
            return GENERAL_ERROR_OUTOFRANGE;
        }
        rfpower = static_cast <uint8_t> (255 * (power - 2.f) / 8.f);
    }
    vector <uint8_t> data { 0x0,0x0 };

    for (size_t i=0; i<2; ++i)
    {
        data[1-i] = 0;
        for (size_t j=0; j<2; ++j)
        {
            uint8_t digit = rfpower % 10;
            switch (j)
            {
            case 0:
                data[1-i] += digit;
                break;
            case 1:
                data[1-i] += digit << 4;
                break;
            }
            rfpower /= 10;
        }
    }

    iretn = usrp_write(handle, 0x14, 0xa, data);
    if (iretn < 0)
    {
        return iretn;
    }

    handle.maxpower = power;
    return 0;
}

int32_t usrp_check_address(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x19, 0x0);
    if (iretn < 0)
    {
        return iretn;
    }

    //	if (handle.address != handle.response[0])
    //	{
    //		return GENERAL_ERROR_OUTOFRANGE;
    //	}

    return 0;
}

int32_t usrp_get_frequency(usrp_handle &handle)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

//    iretn = usrp_write(handle, 0x3);
//    if (iretn < 0)
//    {
//        return iretn;
//    }

//    if (iretn != 5)
//    {
//        return GENERAL_ERROR_OUTOFRANGE;
//    }

//    double frequency = 0.;
//    for (size_t i=0; i<5; ++i)
//    {
//        frequency *= 100.;
//        frequency += 10. * (handle.response[4-i] >> 4) + (handle.response[4-i] % 16);
//    }
//    handle.frequency = frequency;

//    handle.freqband = usrp_freq2band(frequency);
    return iretn;
}

int32_t usrp_get_mode(usrp_handle &handle)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    iretn = usrp_write(handle, 0x4);
    if (iretn < 0)
    {
        return iretn;
    }

    if (iretn != 2)
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }

    switch (handle.response[0])
    {
    case USRP_MODE_AM:
        handle.opmode = DEVICE_RADIO_MODE_AM;
        break;
    case USRP_MODE_FM:
        handle.opmode = DEVICE_RADIO_MODE_FM;
        break;
    case USRP_MODE_LSB:
        handle.opmode = DEVICE_RADIO_MODE_LSB;
        break;
    case USRP_MODE_USB:
        handle.opmode = DEVICE_RADIO_MODE_USB;
        break;
    case USRP_MODE_CW:
    case USRP_MODE_CWR:
        handle.opmode = DEVICE_RADIO_MODE_CW;
        break;
    case USRP_MODE_RTTY:
    case USRP_MODE_RTTYR:
        handle.opmode = DEVICE_RADIO_MODE_RTTY;
        break;
    case USRP_MODE_DV:
        handle.opmode = DEVICE_RADIO_MODE_DV;
        break;
    }
    handle.mode = handle.response[0];
    handle.filtband = handle.response[1];

    iretn = usrp_get_datamode(handle);
    if (iretn < 0)
    {
        return iretn;
    }

    return iretn;
}

int32_t usrp_set_frequency(usrp_handle &handle, double frequency)
{
    int32_t iretn = 0;
    if (iretn < 0)
    {
        return iretn;
    }

    vector <uint8_t> data { 0x0,0x0,0x0,0x0,0x0 };

    if (frequency >= 1e10 || frequency < 0)
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }

    frequency = trunc(frequency);
    for (size_t i=0; i<5; ++i)
    {
        data[i] = 0;
        for (size_t j=0; j<2; ++j)
        {
            uint8_t digit = static_cast <uint8_t>(fmod(frequency, 10.));
            switch (j)
            {
            case 0:
                data[i] += digit;
                break;
            case 1:
                data[i] += digit << 4;
                break;
            }
            frequency = trunc(frequency / 10.);
        }
    }

    iretn = usrp_write(handle, 0x5, data);
    if (iretn < 0)
    {
        return iretn;
    }

    if (frequency < 1.8e6)
    {
        handle.freqband = 14;
    }
    else if (frequency < 2.0e6)
    {
        handle.freqband = 1;
    }
    else if (frequency >= 3.4e6 && frequency < 4.1e6)
    {
        handle.freqband = 2;
    }
    else if (frequency >= 6.9e6 && frequency < 7.5e6)
    {
        handle.freqband = 3;
    }
    else if (frequency >= 9.9e6 && frequency < 10.5e6)
    {
        handle.freqband = 4;
    }
    else if (frequency >= 13.9e6 && frequency < 14.5e6)
    {
        handle.freqband = 5;
    }
    else if (frequency >= 17.9e6 && frequency < 18.5e6)
    {
        handle.freqband = 6;
    }
    else if (frequency >= 20.9e6 && frequency < 21.5e6)
    {
        handle.freqband = 7;
    }
    else if (frequency >= 24.4e6 && frequency < 25.1e6)
    {
        handle.freqband = 8;
    }
    else if (frequency >= 28.0e6 && frequency < 30.0e6)
    {
        handle.freqband = 9;
    }
    else if (frequency >= 50.0e6 && frequency <= 54.0e6)
    {
        handle.freqband = 10;
    }
    else if (frequency >= 108.0e6 && frequency <= 174.0e6)
    {
        handle.freqband = 11;
    }
    else if (frequency >= 420.0e6 && frequency <= 480.0e6)
    {
        handle.freqband = 12;
    }
    else if (frequency >= 1240.0e6 && frequency <1320.0e6)
    {
        handle.freqband = 13;
    }
    else
    {
        handle.freqband = 14;
    }
    handle.frequency = frequency;

    return 0;
}

int32_t usrp_set_mode(usrp_handle &handle, uint8_t opmode)
{
    int32_t iretn;

    iretn = usrp_set_mode(handle, opmode, handle.filtband);
    if (iretn < 0)
    {
        return iretn;
    }

    return 0;

}

int32_t usrp_set_mode(usrp_handle &handle, uint8_t opmode, uint8_t filtband)
{
    int32_t iretn = 0;
    if (iretn < 0)
    {
        return iretn;
    }

    uint8_t datamode = 0;
    uint8_t mode;

    switch (opmode)
    {
    case DEVICE_RADIO_MODE_AMD:
        datamode = 1;
        mode = USRP_MODE_AM;
        break;
    case DEVICE_RADIO_MODE_AM:
        mode = USRP_MODE_AM;
        break;
    case DEVICE_RADIO_MODE_CW:
        mode = USRP_MODE_CW;
        break;
    case DEVICE_RADIO_MODE_CWR:
        mode = USRP_MODE_CWR;
        break;
    case DEVICE_RADIO_MODE_DVD:
        datamode = 1;
        mode = USRP_MODE_DV;
        break;
    case DEVICE_RADIO_MODE_DV:
        mode = USRP_MODE_DV;
        break;
    case DEVICE_RADIO_MODE_FMD:
        datamode = 1;
        mode = USRP_MODE_FM;
        break;
    case DEVICE_RADIO_MODE_FM:
        mode = USRP_MODE_FM;
        break;
    case DEVICE_RADIO_MODE_LSBD:
        datamode = 1;
        mode = USRP_MODE_LSB;
        break;
    case DEVICE_RADIO_MODE_LSB:
        mode = USRP_MODE_LSB;
        break;
    case DEVICE_RADIO_MODE_RTTY:
        mode = USRP_MODE_RTTY;
        break;
    case DEVICE_RADIO_MODE_RTTYR:
        mode = USRP_MODE_RTTYR;
        break;
    case DEVICE_RADIO_MODE_USBD:
        datamode = 1;
        mode = USRP_MODE_USB;
        break;
    case DEVICE_RADIO_MODE_USB:
        mode = USRP_MODE_USB;
        break;
    default:
        return GENERAL_ERROR_OUTOFRANGE;
    }

    if (filtband)
    {
        vector <uint8_t> data { 0x0, 0x0 };
        data[0] = mode;
        data[1] = filtband;
        iretn = usrp_write(handle, 0x6, data);
    }
    else
    {
        vector <uint8_t> data { 0x0 };
        data[0] = mode;
        iretn = usrp_write(handle, 0x6, data);
    }

    if (iretn < 0)
    {
        return iretn;
    }

    handle.filtband = filtband;

    if (handle.datamode != datamode)
    {
        iretn = usrp_set_datamode(handle, datamode);
        if (iretn < 0)
        {
            return iretn;
        }
    }


    handle.opmode = opmode;
    handle.mode = mode;

    return 0;
}

int32_t usrp_set_channel(usrp_handle &handle, uint8_t channelnum)
{
    int32_t iretn = 0;

    switch (channelnum)
    {
    case USRP_CHANNEL_A:
    case USRP_CHANNEL_B:
        {
            iretn = usrp_write(handle, 0x7, 0xd0+channelnum);
        }
        break;
    case USRP_CHANNEL_SWAP:
        {
            iretn = usrp_write(handle, 0x7, 0xb0);
        }
        break;
    default:
        return GENERAL_ERROR_OUTOFRANGE;
    }

    if (iretn < 0)
    {
        return iretn;
    }

    if (channelnum != USRP_CHANNEL_SWAP)
    {
        handle.channelnum = channelnum;
    }

    return 0;
}

int32_t usrp_get_rfgain(usrp_handle &handle)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    iretn = usrp_write(handle, 0x14, 0x2);
    if (iretn != 2)
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }

    handle.rfgain = usrp_byte(handle.response);

    return iretn;
}

int32_t usrp_get_squelch(usrp_handle &handle)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    iretn = usrp_write(handle, 0x14, 0x3);
    if (iretn != 2)
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }
    handle.squelch = usrp_byte(handle.response);

    return 0;
}

int32_t usrp_get_rfpower(usrp_handle &handle)
{
    int32_t iretn = 0;

    if (iretn < 0)
    {
        return iretn;
    }

    iretn = usrp_write(handle, 0x14, 0xa);
    if (iretn != 2)
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }

    float power = usrp_byte(handle.response);

    handle.rfpower = static_cast <uint8_t> (power);
    power /= 255.f;

    if (handle.freqband < 11)
    {
        power = 2.f + power * (handle.mode==USRP_MODE_AM?28.f:98.f);
    }
    else if (handle.freqband < 12)
    {
        power = 2.f + power * 98.f;
    }
    else if (handle.freqband < 13)
    {
        power = 2.f + power * 73.f;
    }
    else if (handle.freqband < 14)
    {
        power = 2.f + power * 8.f;
    }
    handle.maxpower = power;
    return iretn;
}

int32_t usrp_get_smeter(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x15, 0x2);
    if (iretn < 0)
    {
        return iretn;
    }

    float power = usrp_byte(handle.response);

    handle.smeter = static_cast <uint8_t> (power);
    power /= 240.f;

    if (power > 1.f)
    {
        power = 1.;
    }

    if (handle.freqband < 11)
    {
        power = 2.f + power * (handle.mode==USRP_MODE_AM?28.f:98.f);
    }
    else if (handle.freqband < 12)
    {
        power = 2.f + power * 98.f;
    }
    else if (handle.freqband < 13)
    {
        power = 2.f + power * 73.f;
    }
    else if (handle.freqband < 14)
    {
        power = 1.f + power * 8.f;
    }
    handle.powerin = power;

    return 0;
}

int32_t usrp_get_rfmeter(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x15, 0x11);
    if (iretn < 0)
    {
        return iretn;
    }

    float power = usrp_byte(handle.response);

    handle.rfmeter = static_cast <uint8_t> (power);
    if (power <= 141)
    {
        power /= 282.f;
    }
    else
    {
        power /= 215.f;
    }

    if (power > 1.f)
    {
        power = 1.;
    }

    if (handle.freqband < 11)
    {
        power = 2.f + power * (handle.mode==USRP_MODE_AM?28.f:98.f);
    }
    else if (handle.freqband < 12)
    {
        power = 2.f + power * 98.f;
    }
    else if (handle.freqband < 13)
    {
        power = 2.f + power * 73.f;
    }
    else if (handle.freqband < 14)
    {
        power = 1.f + power * 8.f;
    }
    handle.powerout = power;

    return 0;
}

int32_t usrp_get_swrmeter(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x15, 0x12);
    if (iretn < 0)
    {
        return iretn;
    }

    handle.swrmeter = usrp_byte(handle.response);

    return 0;
}

int32_t usrp_get_alcmeter(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x15, 0x13);
    if (iretn < 0)
    {
        return iretn;
    }

    handle.alcmeter = usrp_byte(handle.response);

    return 0;
}

int32_t usrp_get_compmeter(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x15, 0x14);
    if (iretn < 0)
    {
        return iretn;
    }

    usrp_byte(handle.response);

    return 0;
}

int32_t usrp_set_freqband(usrp_handle &handle, uint8_t band)
{
    int32_t iretn = 0;

    if (band > 0 && band < 15)
    {
        vector <uint8_t> data { 0x0, 0x1 };
        if (band < 10)
        {
            data[0] = band;
        }
        else
        {
            data[0] = band + 6;
        }
        iretn = usrp_write(handle, 0x1a, 0x1, data);
    }
    else
    {
        return GENERAL_ERROR_OUTOFRANGE;
    }

    if (iretn < 0)
    {
        return iretn;
    }

    handle.freqband = band;

    return 0;
}

int32_t usrp_get_freqband(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x1a, 0x3);
    if (iretn < 0)
    {
        return iretn;
    }

    if (handle.response[0] < 10)
    {
        handle.bps9600mode = handle.response[0];
    }
    else
    {
        handle.bps9600mode = handle.response[0] - 6;
    }

    return 0;
}

int32_t usrp_set_bps9600mode(usrp_handle &handle, uint8_t mode)
{
    int32_t iretn = 0;

    switch (mode)
    {
    case USRP_9600MODE_OFF:
    case USRP_9600MODE_ON:
        {
            vector <uint8_t> data { 0x0, 0x55, 0x0 };
            data[2] = mode;
            iretn = usrp_write(handle, 0x1a, 0x5, data);
        }
        break;
    default:
        return GENERAL_ERROR_OUTOFRANGE;
    }

    if (iretn < 0)
    {
        return iretn;
    }

    handle.bps9600mode = mode;

    return 0;
}

int32_t usrp_get_bps9600mode(usrp_handle &handle)
{
    int32_t iretn = 0;

    vector <uint8_t> data { 0x0, 0x55 };
    iretn = usrp_write(handle, 0x1a, 0x5, data);
    if (iretn < 0)
    {
        return iretn;
    }

    handle.bps9600mode = handle.response[2];

    return 0;
}

int32_t usrp_set_datamode(usrp_handle &handle, uint8_t mode)
{
    int32_t iretn = 0;

    switch (mode)
    {
    case USRP_DATAMODE_ON:
    case USRP_DATAMODE_OFF:
        {
            vector <uint8_t> data { 0x0, 0x0 };
            data[0] = mode;
            if (mode == USRP_DATAMODE_ON)
            {
                data[1] = handle.filtband;
                if (data[1] == 0 || data[1] > 3)
                {
                    data[1] = 1;
                }
            }
            iretn = usrp_write(handle, 0x1a, 0x6, data);
            if (iretn < 0)
            {
                return iretn;
            }
            iretn = usrp_set_bps9600mode(handle, USRP_9600MODE_ON);
        }
        break;
    default:
        return GENERAL_ERROR_OUTOFRANGE;
    }

    if (iretn < 0)
    {
        return iretn;
    }

    handle.datamode = mode;
    handle.opmode = static_cast <uint8_t>((handle.opmode >> 1) << 1) + mode;

    return 0;
}

int32_t usrp_get_datamode(usrp_handle &handle)
{
    int32_t iretn = 0;

    iretn = usrp_write(handle, 0x1a, 0x6);
    if (iretn < 0)
    {
        return iretn;
    }

    handle.datamode = handle.response[0];
    handle.opmode = static_cast <uint8_t>((handle.opmode >> 1) << 1) + handle.response[0];

    return 0;
}