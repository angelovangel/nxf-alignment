#!/usr/bin/env python3
import sys
import re

mapping = {
    'C+m': '5mC',
    'C+h': '5hmC',
    'C+f': '5fC',
    'C+c': '5caC',
    'A+a': '6mA',
    '+b': 'Abasic',
    'C+e': '5eC'
}

def main():
    detected = set()
    # Read all input and split by whitespace and semicolons
    input_data = sys.stdin.read()
    normalized = input_data.replace(';', ' ')
    tags = normalized.split()
    
    for tag in tags:
        tag = tag.strip()
        if not tag:
            continue
            
        found = False
        for pattern, name in mapping.items():
            if pattern in tag: 
                detected.add(name)
                found = True
        
        if not found:
            pass
    
    if not detected:
        sys.stdout.write("-")
    else:
        # Concatenate found mods
        sys.stdout.write(", ".join(sorted(list(detected))))

if __name__ == "__main__":
    main()
