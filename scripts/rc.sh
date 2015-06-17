#!/bin/sh

TEMP="/tmp/.tmp_$(date +%s)"

cut -f 1 $1 | sed 's/^/>\n/' > $TEMP
fastools reverse $TEMP - | grep -v '>'

rm $TEMP
