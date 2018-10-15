from abc import abstractmethod, ABC


class AnsatzGenerator(ABC):

    @abstractmethod
    def generate(self, inputParams):
        pass
