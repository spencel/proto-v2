#!/usr/bin/python

# Data structure of workflow config files
# {
#   <label>: {
#       "command": <command>,
#       "arguments": {
#           <cli-option>: <value>,
#           ...
#       },
#   ...
#   }
# }
# Example
# ...
# "Process PCR Diagnostic Tests": {
#     "command": "process_diagnostic_tests",
#     "arguments": {
#         "--dirpath": "./data/mock/diagnostic-kit-definition-files/pcr",
#         "--test-type": "PCR"
#         }
#     },
# "Process LAMP Diagnostic Tests": {
#     "command": "process_diagnostic_tests",
#     "arguments": {
#         "--dirpath": "./data/mock/diagnostic-kit-definition-files/lamp",
#         "--test-type": "LAMP"
#     }
#  },
#  ...

import json
import logging

import commands
import config

class Workflow():
    

    def __init__(
        self,
        config,
        first_step=1,
        skip_steps=[],
        last_step=999999999,
        skip_commands=[""],
        skip_labels=[""]
    ):
        # Handle config arg
        # Case 1: Arg is a filepath
        if isinstance(config, str):
            self.config = json.load(open(config))
        # Case 2: Assume arg is already a dictionary
        else:
            self.config = config

        # Handle other args
        self.first_step = first_step
        self.skip_steps = skip_steps
        self.last_step = last_step
        self.skip_commands = skip_commands
        self.skip_labels = skip_labels

    def exec(self):
        scope_pref = f"{Workflow.__name__}.{Workflow.exec.__name__}()"
        logging.info(f"{scope_pref}: started")

        # Configs
        # Note used because commands are not executed from CLI
        python_exe = config.paths.python_exe

        # Configs: alias properties
        wf_config = self.config
        current_step = self.first_step
        skip_steps = self.skip_steps
        last_step = self.last_step
        skip_commands = self.skip_commands
        skip_labels = self.skip_labels
        # logging.debug(f"{scope_pref}: skip_labels: {skip_labels}")

        # Exit with error if config contains duplicate keys, ie, labels
        labels = wf_config.keys()
        labels_set = set(labels)
        if len(labels) > len(labels_set):
            raise Exception("ERROR: Workflow config json file contains duplicate labels.")

        # Set step
        step_counter = 1
        for label in wf_config:
            command_name = wf_config[label]["command"]
            # logging.debug(f"{scope_pref}: command_name: {command_name}")
            # Create command statement
            # Get command
            str_command = f"commands.{command_name}(\n"
            # Get arguments
            args = dict()
            if "arguments" in wf_config[label]:
                args = wf_config[label]["arguments"]
            # Handles args
            for arg in args:
                # Get argument value
                arg_value = args[arg]
                # Give argument value quotes if it's a string data type
                if isinstance(arg_value, str):
                    arg_value = f"\"{arg_value}\""
                # Convert arg name to function parameter name if necessary
                arg = arg.lstrip("-").replace("-", "_")
                str_command = f"{str_command}\t{arg}={arg_value},"
            # Remove trailing "," and close parenthesis
            str_command = f"{str_command.rstrip(',')}\n)"
            # logging.debug(f"{scope_pref}: str_command B: {str_command}")

            if (
                current_step == step_counter
                and current_step not in skip_steps
                and command_name not in skip_commands
                and label not in skip_labels
                ):
                logging.info(f"{scope_pref}: starting step {step_counter}: {label}")
                eval(str_command)
                current_step += 1
            else:
                logging.info(f"{scope_pref}: skipped step {step_counter}: {label}")

            if (
                current_step <= step_counter
                and (
                    current_step in skip_steps
                    or command_name in skip_commands
                    or label in skip_labels
                )
                ):
                current_step += 1

            if current_step > last_step:
                logging.info(f"{scope_pref}: stopped after step {step_counter}: {label}")
                break

            step_counter += 1
            # logging.debug(f"{scope_pref}: current_step, step_counter: {current_step}, {step_counter}")

        logging.info(f"{scope_pref}: end of workflow")
