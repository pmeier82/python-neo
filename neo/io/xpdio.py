# encoding: utf-8
"""
xpdio
======
.xpd binary files from Logothetis lab

Classes
-------

_XPD_TRIAL_HEADER   - class to read a trial header
_XPD_CHANNEL_HEADER - class to read a channel header
XpdIO               - class to read one .xpd file
"""

##---IMPORTS

import os
import scipy as sp
import quantities as pq
import warnings
from struct import Struct
from .baseio import BaseIO
from .tools import create_many_to_one_relationship
from ..core import Block, Segment, AnalogSignalArray, RecordingChannelGroup

##---CONSTANTS

_CONV_I = Struct('I')
_CONV_H = Struct('H')
_CONV_B = Struct('B')
_CONV_d = Struct('d')
_CONV_s = Struct('s')
_CONV_7s = Struct('7s')
_CONV_255s = Struct('255s')

##---CLASSES

class _XPD_TRIAL_HEADER(object):
    """XPD trial header

    ===== =============================
    Size  XPD Trial Header Member
    ===== =============================
    H     header size
    I     data size
    7s    name of generating program
    B     version identifier = 86 / 'V'
    H     software version
    B     high ver. = 0
    B     low ver. = 12
    I     trial no.
    B     stimulus code
    B     error byte
    4H    time stamp
    255s  comment
    255s  additional comment
    Total 542 bytes
    ===== =============================
    """

    SIZE = 542

    def __init__(self, fp):
        """
        :type fp: FileIO
        :param fp: open file at seek(header_start)
        """

        # read buffer
        buf = fp.read(self.SIZE)

        # extract trial header information
        self.header_size = _CONV_H.unpack(buf[0:2])[0]
        self.data_size = _CONV_I.unpack(buf[2:6])[0]
        self.name = _CONV_7s.unpack(buf[6:13])[0]
        self.version = _CONV_s.unpack(buf[13:14])[0]
        self.sw_build_no = _CONV_H.unpack(buf[14:16])[0]
        self.high_ver = _CONV_B.unpack(buf[16:17])[0]
        self.low_ver = _CONV_B.unpack(buf[17:18])[0]
        self.trial_no = _CONV_I.unpack(buf[18:22])[0]
        self.stm_code = _CONV_B.unpack(buf[22:23])[0]
        self.err_code = _CONV_B.unpack(buf[23:24])[0]
        self.timestamp = (_CONV_H.unpack(buf[24:26])[0],
                          _CONV_H.unpack(buf[26:28])[0],
                          _CONV_H.unpack(buf[28:30])[0],
                          _CONV_H.unpack(buf[30:32])[0])
        self.comment = _CONV_255s.unpack(buf[32:287])[0]
        self.comment = self.comment[:self.comment.find('\x00')]
        self.add_comment = _CONV_255s.unpack(buf[287:542])[0]
        self.add_comment = self.add_comment[:self.add_comment.find('\x00')]

    def __str__(self):
        rval = self.__repr__()
        rval += '\nheader_size\t\t%d\n' % self.header_size
        rval += 'data_size\t\t%d\n' % self.data_size
        rval += 'name\t\t\t%s\n' % self.name
        rval += 'sw_ver\t\t\t%s%s\n' % (self.version, self.sw_build_no)
        rval += 'pt_ver\t\t(%d, %d)\n' % (self.high_ver, self.low_ver)
        rval += 'trial_no\t\t%d\n' % self.trial_no
        rval += 'stm_code\t\t%d\n' % self.stm_code
        rval += 'err_code\t\t%d\n' % self.err_code
        rval += 'timestamp\t\t%s\n' % str(self.timestamp)[1:-1]
        rval += 'comment\t\t\t%s\n' % self.comment
        rval += 'add_comment\t\t%s\n' % self.add_comment
        return rval

    def __unicode__(self):
        return unicode(self.__str__())


class _XPD_CHANNEL_HEADER(object):
    """XPD channel header

        ===== =====================
        Size  XPD Channel Header
        ===== =====================
        I     channel no.
        d     sample rate in kHz
        d     offset-x of channel
        I     data length in samples
        Total 24 bytes
        ===== =====================
    """

    SIZE = 24

    def __init__(self, fp):
        """
        :type fp: FileIO
        :param fp: open file at seek(header_start)
        """

        # read buffer
        buf = fp.read(self.SIZE)

        # extract channel header information
        self.channel_no = _CONV_I.unpack(buf[0:4])[0]
        self.sample_rate = _CONV_d.unpack(buf[4:12])[0]
        self.x_offset = _CONV_d.unpack(buf[12:20])[0]
        self.n_sample = _CONV_I.unpack(buf[20:24])[0]
        self.data_offset = fp.tell()

    def __str__(self):
        rval = self.__repr__()
        rval += '\nchannel_no\t\t%d\n' % self.channel_no
        rval += 'sample_rate\t\t%f\n' % self.sample_rate
        rval += 'x_offset\t\t%f\n' % self.x_offset
        rval += 'n_sample\t\t%d\n' % self.n_sample
        rval += 'data_offset\t\t%d' % self.data_offset
        return rval

    def __unicode__(self):
        return unicode(self.__str__())


class XpdIO(BaseIO):
    """yadda yadda"""

    ## class variables

    is_readable = True
    is_writable = False
    is_streameable = False
    has_header = True

    supported_objects = [Segment]
    readable_objects = [Segment]
    writeable_objects = []

    ## gui info

    read_params = None
    write_params = None
    name = 'XpdIO'
    description = 'IO for the .xpd format from Logothetis lab, MPI Tuebingen'
    extensions = ['xpd']
    mode = 'file'

    ## constructor

    def __init__(self, filename=None, **kargs):
        """"""

        # super
        super(XpdIO, self).__init__(filename)

        # members
        self._file = None
        self._header = None
        self._ch_analog = {}
        self._ch_digital = {}
        self._ch_event = {}

        # open file handle
        self._open()

    def __del__(self):
        self._close()
        del self._ch_analog, self._ch_digital, self._ch_event

    ## private interface

    def _open(self):
        """opens the file handle"""
        # open file
        self._file = open(self.filename, 'rb')

        # trial header
        if _CONV_H.unpack(self._file.read(2))[0] != 120:
            self._close()
            raise IOError('unexpected input while reading trial announce for file %s' % self._file.name)
        self._header = _XPD_TRIAL_HEADER(self._file)

        # analog channel headers
        if _CONV_H.unpack(self._file.read(2))[0] != 123:
            self._close()
            raise IOError('unexpected input while reading analog channel announce for file %s' % self._file.name)
        n_ch_a = _CONV_I.unpack(self._file.read(4))[0]
        self._ch_analog = {}
        for _ in xrange(n_ch_a):
            ach = _XPD_CHANNEL_HEADER(self._file)
            self._ch_analog[ach.channel_no] = ach
            self._file.seek(ach.n_sample * 2, 1)

        # digital channel headers
        if _CONV_H.unpack(self._file.read(2))[0] != 121:
            self._close()
            raise IOError('unexpected input while reading digital channel announce for file %s' % self._file.name)
        n_ch_d = _CONV_I.unpack(self._file.read(4))[0]
        self._ch_digital = {}
        for _ in xrange(n_ch_d):
            dch = _XPD_CHANNEL_HEADER(self._file)
            self._ch_digital[dch.channel_no] = dch
            self._file.seek(dch.n_sample * 4, 1)

        # event channel headers
        if _CONV_H.unpack(self._file.read(2))[0] != 122:
            self._close()
            raise IOError('unexpected input while reading event channel header for file %s' % self._file.name)
        n_ch_e = _CONV_I.unpack(self._file.read(4))[0]
        self._ch_event = {}
        for _ in xrange(n_ch_e):
            ech = _XPD_CHANNEL_HEADER(self._file)
            self._ch_event[ech.channel_no] = ech
            self._file.seek(ech.n_sample * 4, 1)

        # rewind
        self._file.seek(0)

    def _close(self):
        """close the file handle and stop operations"""
        if self._file:
            if not self._file.closed:
                self._file.close()

    def _get_channels(self, chans, kind='a'):
        """returns a numpy array of the channel data"""

        # init and check
        if kind not in ['a', 'd', 'e']:
            raise ValueError('unknown kind (%s), accepts: \'a\', \'d\' and \'e\'!')
        ch = {'a': self._ch_analog,
              'd': self._ch_digital,
              'e': self._ch_event}[kind]
        dtype = {'a': sp.int16,
                 'd': sp.int32,
                 'e': sp.int32}[kind]
        if not chans:
            raise ValueError('chans not valid!')
        n_sample = 0
        for c in chans:
            try:
                if ch[c].n_sample > n_sample:
                    n_sample = ch[c].n_sample
            except KeyError:
                warnings.warn('invalid channel: %s' % c, RuntimeWarning)
        if n_sample == 0:
            warnings.warn('no data for chans %s' % str(chans), RuntimeWarning)
            return sp.zeros((n_sample, len(chans)), dtype=dtype)

        # collect data
        rval = sp.zeros((n_sample, len(chans)), dtype=dtype)
        for i, c in enumerate(chans):
            # load from file
            try:
                header = ch[c]
                self._file.seek(header.data_offset)
                load_item = sp.fromfile(file=self._file, count=header.n_sample, dtype=dtype)
            except KeyError:
                warnings.warn('invalid channel: %s!' % c, RuntimeWarning)
                load_item = sp.zeros(n_sample, dtype=dtype)

            # unfortunately the channel shapes are not always consistent
            # across the tetrode. but we preallocate space such that the
            # largest item may fit in the buffer.
            rval[:load_item.size, i] = load_item

        # return stuff
        return rval

    def _available_tetrodes(self):
        """yields the set of available tetrodes"""

        return [k for k in self._ch_analog.keys() if 1 <= k <= 16]

        ## neo components read interface

    def read(self, **kwargs):
        return self.read_block(**kwargs)

    def read_block(self, **kargs):
        """read the block

        :keywords: none
        """

        # init
        trial_name = os.path.basename(self._file.name).split('.')[-2]
        exp_name = trial_name[:4]
        trial_name = trial_name[4:]

        # create block
        blk = Block(name=exp_name)

        # create segment
        seg = Segment(
            name=trial_name,
            file_origin=self._file.name)
        # TODO: annotations from trial header
        blk.segments.append(seg)

        # care for tetrodes
        for t in self._available_tetrodes():
            tet = RecordingChannelGroup(name='Tet::%d' % t)
            blk.recordingchannelgroups.append(tet)
            # TODO: annotations from trial header or xpt?
            tet_cs = [t + i * 16 for i in xrange(4)]
            sr = None
            for c in tet_cs:
                try:
                    sr = self._ch_analog[c].sample_rate
                    break
                except KeyError:
                    continue
            if not sr:
                warnings.warn('no sample rate found! recording group may be invalid!', RuntimeWarning)
                sr = 32.0
            tet_data_raw = self._get_channels(tet_cs, kind='a')
            tet_data = AnalogSignalArray(
                signal=tet_data_raw,
                sampling_rate=sr * pq.kHz,
                units='uV')
            tet.analogsignalarrays.append(tet_data)

        # return
        create_many_to_one_relationship(blk)
        return blk

##---MAIN

if __name__ == '__main__':
    pass
