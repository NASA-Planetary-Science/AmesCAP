#!/usr/bin/env python3
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        description='Welcome to AMESCAP!',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Add arguments
    parser.add_argument('command', nargs='?', default='help',
                       help='Command to execute (use "help" for more information)')

    args = parser.parse_args()

    # If no arguments or help is requested
    if args.command == 'help' or '-h' in sys.argv or '--help' in sys.argv:
        print("""
Welcome to AMESCAP!
------------------

Available commands:
  MarsCalendar.py  - [Brief description]
  MarsFiles.py     - [Brief description]
  MarsFormat.py    - [Brief description]
  MarsInterp.py    - [Brief description]
  MarsPlot.py      - [Brief description]
  MarsPull.py      - [Brief description]
  MarsVars.py      - [Brief description]

For detailed help on each command, use:
  <command> -h
  Example: MarsVars.py -h

[Add any additional information about AMESCAP here]
        """)
        return 0

if __name__ == '__main__':
    sys.exit(main())