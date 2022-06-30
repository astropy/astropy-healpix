import astropy_healpix
from abc import abstractmethod
from astropy.coordinates.representation import (
    REPRESENTATION_CLASSES,
    DUPLICATE_REPRESENTATIONS,
    BaseRepresentation,
    CartesianRepresentation,
)
from astropy.coordinates import Distance
import astropy.units as u


class HEALPixRepresentationBase(BaseRepresentation):

    attr_classes = {
        "indices": u.Quantity,
        "dx": u.Quantity,
        "dy": u.Quantity,
    }

    _order: int
    _nside: int
    _dx: u.Quantity
    _dy: u.Quantity

    @classmethod
    def get_subclass(cls, *, order, nside):
        name: str = f"HEALPix{order.capitalize()}{nside}"
        if "Unit" in cls.__name__:
            name += "Unit"

        if name.lower() in REPRESENTATION_CLASSES:
            return REPRESENTATION_CLASSES[name.lower()]
        elif "abc." + name in REPRESENTATION_CLASSES:
            return REPRESENTATION_CLASSES["abc." + name]

        return type(name + "Representation", (cls,), {}, order=order, nside=nside)

    @classmethod
    @property
    @abstractmethod
    def order(cls):
        return cls._order

    @classmethod
    @property
    @abstractmethod
    def nside(cls):
        return cls._nside

    def __init__(self, *args, differentials=None, copy=True):
        super().__init__(*args, differentials=differentials, copy=copy)
        self._indices = (self._indices << u.one).astype(int)
        self._dx <<= u.one
        self._dy <<= u.one

    @property
    def indices(self):
        """HEALPix indices."""
        return self._indices

    @property
    def dx(self):
        """HEALPix x offset."""
        return self._dx

    @property
    def dy(self):
        """HEALPix y offset."""
        return self._dy


class HEALPixBaseRepresentation(HEALPixRepresentationBase):

    attr_classes = {
        "indices": u.Quantity,
        "distance": Distance,
        "dx": u.Quantity,
        "dy": u.Quantity,
    }

    _distance: Distance

    def __init_subclass__(cls, *, order: str, nside: int, **kwargs):
        super().__init_subclass__(**kwargs)

        astropy_healpix.core._validate_nside(nside)
        cls._nside = int(nside)

        astropy_healpix.core._validate_order(order)
        cls._order = order

    @classmethod
    @property
    @abstractmethod
    def order(cls):
        return cls._order

    @classmethod
    @property
    @abstractmethod
    def nside(cls):
        return cls._nside

    def __init__(
        self, indices, distance, dx=0.5, dy=0.5, differentials=None, copy=True
    ):
        super().__init__(
            indices, distance, dx, dy, differentials=differentials, copy=copy
        )

    @property
    def distance(self):
        """The distance from the origin to the point(s)."""
        return self._distance

    def to_cartesian(self):
        # TODO! register astropy_healpix.healpix_to_xyz with Quantity
        xyz = astropy_healpix.healpix_to_xyz(
            self.indices.value,
            self.nside,
            dx=self.dx.value,
            dy=self.dy.value,
            order=self.order,
        )
        return CartesianRepresentation(*xyz * self.distance)

    @classmethod
    def from_cartesian(cls, cart):
        distance = cart.norm()
        indices, dx, dy = astropy_healpix.xyz_to_healpix(
            *cart.xyz.value, nside=cls.nside, return_offsets=True, order=cls.order
        )
        return cls(indices, distance, dx=dx, dy=dy)


class HEALPixBaseUnitRepresentation(HEALPixRepresentationBase):
    def __init_subclass__(cls, *, order: str, nside: int, **kwargs):
        super().__init_subclass__(**kwargs)

        astropy_healpix.core._validate_nside(nside)
        cls._nside = int(nside)

        astropy_healpix.core._validate_order(order)
        cls._order = order

    def __init__(self, indices, dx=0.5, dy=0.5, differentials=None, copy=True):
        super().__init__(indices, dx, dy, differentials=differentials, copy=copy)

    @classmethod
    @property
    @abstractmethod
    def order(cls):
        return cls._order

    @classmethod
    @property
    @abstractmethod
    def nside(cls):
        return cls._nside

    def to_cartesian(self):
        xyz = astropy_healpix.healpix_to_xyz(
            self.indices.value,
            self.nside,
            dx=self.dx.value,
            dy=self.dy.value,
            order=self.order,
        )
        return CartesianRepresentation(*xyz)

    @classmethod
    def from_cartesian(cls, cart):
        distance = cart.norm()
        indices, dx, dy = astropy_healpix.xyz_to_healpix(
            *cart.xyz.value, nside=cls.nside, return_offsets=True, order=cls.order
        )
        return cls(indices, dx=dx, dy=dy)
