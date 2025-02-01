import torch

from args_from_cli import ArgsParams
from models import Encoder as Encoder


class PartialDeformationEncoderBase(torch.nn.Module):
    """
    Abstract class with API for encoders
    """

    def __init__(self):
        super(PartialDeformationEncoderBase, self).__init__()

    def to_load(self):
        """
        :return: string, name of field to be loaded for the batch, e.g., if the returned value is "pie" then the data loader
        wiil attempt to load the file "pie.npy" for each example in the batch. This field will be then fed as an argument
        to the to_code function
        """
        raise NotImplementedError

    def encode(self, loaded):
        """
        Return the code for the given example (represented by "loaded")
        :param loaded: the content of the loaded field requested in to_load()
        :return: the code for the loaded example
        """
        raise NotImplementedError


class _LoadedDataEncoder(PartialDeformationEncoderBase):
    """
    Loads a given key_name as an npy file, and retuns it stacked as a 1D vector
    """

    def __init__(self, key_name):
        super(_LoadedDataEncoder, self).__init__()
        self.__key_name = key_name

    def to_load(self):
        return self.__key_name

    def encode(self, loaded):
        return loaded


class _PointNetEncoder(PartialDeformationEncoderBase):
    """
    encodes the "samples.npy" field via a pointnet
    """

    def __init__(self, output_dim, normalization):
        super(_PointNetEncoder, self).__init__()
        self.__pointnet = Encoder(output_dim, self.__get_in_dim(), normalization=normalization)

    def __get_in_dim(self):
        in_dim = 3
        in_dim += 3
        return in_dim

    def to_load(self):
        load_list = ["samples"]
        load_list.append("samples_normals")
        return load_list

    def encode(self, loaded):
        r = torch.cat(loaded, dim=2)
        r = self.__pointnet(r)
        return r

    def encode_batch(self, loaded):
        r_list = []
        for i in range(len(loaded)):
            r_list.append(torch.cat(loaded[i], dim=2))
        r = torch.cat(r_list, dim=0)
        r = self.__pointnet(r)
        return r


class DeformationEncoder(torch.nn.Module):
    """
    Class for encoding deformations
    """

    SOURCE = True
    TARGET = False

    def __init__(self, args: ArgsParams):
        super(DeformationEncoder, self).__init__()
        self.__keys_to_load = {
            DeformationEncoder.SOURCE: {},
            DeformationEncoder.TARGET: {},
        }
        self.__generators = {
            DeformationEncoder.SOURCE: [],
            DeformationEncoder.TARGET: [],
        }
        self.__module_list = torch.nn.ModuleList()
        self.__ttype = None
        self.__args = args

    def __add_generator(self, source_or_target: bool, generator: PartialDeformationEncoderBase):
        """
        add the given encoder
        :param source_or_target: add it to encode the source, or the target?
        :param generator: the encoder
        """
        keys = generator.to_load()
        if isinstance(keys, str):
            keys = [keys]
        for key in keys:
            self.__keys_to_load[source_or_target][key] = True
        self.__generators[source_or_target].append(generator)
        self.__module_list.append(generator)

    def add_pointnet(self, code_length: int, source: bool, target: bool):
        """
        Add a pointnet encoder
        :param code_length: the desired code length of PN's output
        :param source: true/false -- apply the PN to the source mesh?
        :param target: true/false -- apply the PN to the target mesh?
        """
        encoder = _PointNetEncoder(
            code_length,
            normalization=self.__args.pointnet_layer_normalization,
        )
        if self.__ttype is not None:
            encoder.type(self.type())
        to_add = []
        if source:
            to_add.append(True)
        if target:
            to_add.append(False)
        assert len(to_add) > 0
        for val in to_add:
            self.__add_generator(val, encoder)

    def add_loader(self, source_or_target: bool, field_name):
        """
        Add an encoder that loads the given field name ("dog" will load "dog.npy" for the given example in the batch)[
        :param source_or_target: this loader is applied to source or target?
        :param field_name: the fieldname to load
        :return:
        """
        gem = _LoadedDataEncoder(field_name)
        self.__add_generator(source_or_target, gem)

    def get_keys_to_load(self, source_or_target: bool):
        """
        get all keys that are needed to load from disk
        :param source_or_target: give keys for source, or target
        :return: list of strings for all keys
        """
        return self.__keys_to_load[source_or_target]

    def __get_partial_code(self, the_obj, source_or_target: bool):
        """
        given a batch object, process and return the partial code representing the source/target
        :param batch: Batch object, for which codes are to be computed
        :param source_or_target: return partial code for source, or target?
        :return:
        """
        codes = []
        for generator in self.__generators[source_or_target]:
            keys = generator.to_load()
            if isinstance(keys, str):
                code = the_obj.get_loaded_data(keys)
            else:
                assert isinstance(keys, list)
                code = []
                for key in keys:
                    code.append(the_obj.get_loaded_data(key))
            code = generator.encode(code)
            code = code.view(code.shape[0], -1)
            codes.append(code)
        if len(codes) == 0:
            return None
        ret = torch.cat(codes, dim=1)
        assert ret.shape[0] == codes[0].shape[0]
        return ret

    def __get_partial_codes_batch(self, list_objs, source_or_target: bool):
        """
        given a batch object, process and return the partial code representing the source/target
        :param batch: Batch object, for which codes are to be computed
        :param source_or_target: return partial code for source, or target?
        :return:
        """
        codes = []
        list_code = []
        for generator in self.__generators[source_or_target]:
            for obj in list_objs:
                keys = generator.to_load()
                if isinstance(keys, str):
                    code = obj.get_loaded_data(keys)
                else:
                    assert isinstance(keys, list)
                    code = []
                    for key in keys:
                        code.append(obj.get_loaded_data(key))
                list_code.append(code)
            list_code = generator.encode_batch(list_code)
            list_code = list_code.view(list_code.shape[0], -1)
            codes.append(list_code)

        if len(codes) == 0:
            return None
        ret = torch.cat(codes, dim=1)
        assert ret.shape[0] == codes[0].shape[0]
        return ret

    def encode_deformation(self, source):
        """
        get the code for the deformations in a batch
        :param batch: Batch object
        :return: a batch of codes representing the deformations in the batch
        """
        s = self.__get_partial_codes_batch(source, DeformationEncoder.SOURCE)
        return s

    def get_code_length(self, dataset):
        """
        return the length of the code for a given example. Since this is decided at run time, this loads one batch from
        the dataset and checks the code length direclty
        :param dataset: Dataset object that returns batches
        :return: integer specifying code lenth
        """
        b1 = dataset[0]

        c = self.encode_deformation((b1[0],))
        return c.shape[1]
