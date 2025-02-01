import os
import warnings
from abc import ABC, abstractmethod

import args_from_cli
from DeformationEncoder import DeformationEncoder

from train_loop import main


class Experiment(ABC):
    """
    base class for experiments
    """

    def __init__(self, name, description):
        self.net = None
        self.name = name
        self.description = description

    def modify_args(self, args):
        """
        called before setting args, to enable modifying args
        :param args: the original args
        :return: the args, with any modification that should be included
        """
        return args

    def get_encoder(self, args):
        """
        initialize the encoder for this experiment and return it
        :param args: the cli args
        :return: the encoder, initialize
        """
        args = self.modify_args(args)
        encoder = DeformationEncoder(args)
        self.init_encoder(encoder, args)
        return encoder

    @abstractmethod
    def init_encoder(self, encoder, args):
        """
        abstract method that should be overridden to init the encoder object
        :return: DeformationEncoder object
        """
        pass

    def get_args_and_train(self):
        if self.net is not None:
            warnings.warn(
                "seems like you loaded a network, but are now running training -- FYI, the loaded network is not being used in training (you need to specify the checkpoint in CLI"
            )
        args = args_from_cli.parse_args()

        args = self.modify_args(args)
        self.args = args
        print(f"starting training with args: {args}")

        if self.args.mode == "test":
            if self.args.no_rotational_alignment:
                gen = os.path.join(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "network weights",
                    "network_unaligned_objects.ckpt",
                )
            else:
                gen = os.path.join(
                    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "network weights",
                    "network_aligned_objects.ckpt",
                )
            print(f"************************** STARTING TEST USING CHECKPOINT {gen}")
        else:
            gen = self.get_encoder(args)

        main(gen, self.name, args)
