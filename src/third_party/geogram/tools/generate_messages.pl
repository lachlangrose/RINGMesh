#!/usr/bin/perl
#
# This script generates encrypted messages printed at Vorpaline startup.
# It generates the message file in the Vorpaline source tree directly
#
use strict;
use warnings FATAL => 'all';
use POSIX();
use FindBin();

#! Message file
my $MESSAGE_FILE = "$FindBin::Bin/../src/lib/vorpalib/license/messages.h";

#! Current year
my $YEAR = POSIX::strftime('%Y', localtime());

#!
# \brief Message table
# \details Each message is specified by 3 or 4 fields:
# -# msg_name - the symbolic name of the message
# -# msg - the text of message to encode
# -# key_name - the name of the encryption key (or its value)
# -# key_offset - optional index of the first key char (default = 0)
#
my $MESSAGES = [
    ["CR_1", "\n", "VOR_MESSAGE_KEY", 9],
    ["CR_2", "\n", "\x00"],
    ["BANNER_VORPALINE", "\\V (O |R |P /A |L |I |N |E", "VOR_MESSAGE_KEY"],
    ["BANNER_COPYRIGHT", "(C) Bruno Levy, INRIA - ALICE, 2012-$YEAR", "VOR_MESSAGE_KEY"],
    ["COMPILED", "Compiled: ", "VOR_MESSAGE_KEY"],
    ["EVALUATION_LICENSE_FOR", "Evaluation license for: ", "VOR_MESSAGE_KEY", 4],
    ["LICENSE_IS", "License is ", "VOR_MESSAGE_KEY"],
    ["PERSONAL_AND_NON_TRANSFERABLE", "personal and non-transferable", "VOR_MESSAGE_KEY"],
    ["EXPIRED_LICENSE", "Expired license", "VOR_MESSAGE_KEY"],
    ["EXPIRES_IN", "Expires in ", "VOR_MESSAGE_KEY"],
    ["DAYS", " days\n", "VOR_MESSAGE_KEY"],
    ["EXPIRES_TODAY", "Expires today.", "VOR_MESSAGE_KEY", 3],

    ["HOME", "HOME", "VOR_FILE_KEY"],

    ["PROC_1", "/proc/1", "VOR_FILE_KEY"],
    ["PROC_SELF", "/proc/self", "VOR_FILE_KEY"],
    ["TMP", "/tmp", "VOR_FILE_KEY"],
    ["VAR_LOG", "/var/log", "VOR_FILE_KEY"],
];

#!
# \brief Encryption keys
# \details TODO: instead of duplicating the key values,
# we should read the keys from license.h directly 
#
my %KEYS = (
    "VOR_FILE_KEY"    => "ATUH\0VUUUUUUUH\0T\$ H+C\0",
    "VOR_MESSAGE_KEY" => "AUATU\0T\$.H\0\\\$\@A\0VUUUUUUUH\0T\$ H\0H+C\0ATUH\0",
);

generate_messages();
exit(0);

##############################################################################
#!
# \brief Generates the text of the message file
#
sub generate_messages {

    my $text = <<END;
/*
 * Encrypted messages for license usage
 * This file is automatically generated from tools/generate_messages.pl
 * Do not edit manually
 */

#ifndef __VORPALINE_LICENSE_MESSAGES__
#define __VORPALINE_LICENSE_MESSAGES__
END

    foreach my $message (@$MESSAGES) {
        $text .= encode_message(@$message);
    }

    $text .= "\n#endif\n\n";

    save_file($MESSAGE_FILE, $text);
}

##############################################################################
#!
# \brief Encodes a single message
# \param[in] msg_name the symbolic name of the message
# \param[in] msg the text of message to encode
# \param[in] key_name the name of the encryption key (or its value)
# \param[in] key_offset optional index of the first key char (default = 0)
#
sub encode_message {
    my($msg_name, $msg, $key_name, $key_offset) = @_;

    #print STDERR "*** Encoding: $msg_name, $msg, $key_name, ", $key_offset || 0, "\n";

    # Get the key value in the key table

    my $key = $KEYS{$key_name};
    if( not defined $key ) {
        # The key name is not in the table, so we take it as a literal key value
        $key = $key_name;

        # Encode the key name has a hexadecimal string
	$key_name =~ s{.}{ sprintf("\\x%x", ord($&)) }ges;
	$key_name = "\"$key_name\"";
    }

    # Check if a the key offset is defined

    if( defined($key_offset) ) {
        $key_name .= "+$key_offset";
    }

    my $encoded_msg = encode_string($msg, $key, $key_offset);

    return <<END

/* Message: '$msg' */
#define VOR_MSG_${msg_name}_KEY $key_name
#define VOR_MSG_${msg_name} \"$encoded_msg\"
END
}

##############################################################################
#!
# \brief Encode the message by XORing message chars with key chars
# \param[in] msg the message to encode
# \param[in] key the encryption key
# \param[in] key_index index of the first key char (default = 0)
#
sub encode_string {
    my($msg, $key, $key_index) = @_;

    if( not defined($key_index) ) {
        $key_index = 0;
    }

    my $msg_length = length($msg);
    my @key_chars = split(//, $key);

    $msg =~ s{.}{ sprintf("\\x%x", ord($&) ^ ord($key_chars[$key_index++])) }ges;
    return sprintf("\\x%x%s", $msg_length, $msg);
}

##############################################################################
#!
# \brief Writes a text to a file
# \param[in] file path to the file to save
# \param[in] text text to write to \p file
#
sub save_file {
    my($file, $text) = @_;

    print "Generating file $file\n";

    if( not open(FILE, '>', $file) ) {
        die "Error: failed to save file $file: $!\n";
    }

    print FILE $text;
    close(FILE);
}



