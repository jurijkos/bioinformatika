#!/usr/bin/env bash
touch "$3"
echo ">gene from $1" >> "$3"
cat "$1" >> "$3"
echo ">gene from $2" >> "$3"
cat "$2" >> "$3"
