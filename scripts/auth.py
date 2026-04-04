"""
Authentication Module for AlphaGWAS API.

Provides:
- JWT token-based authentication
- API key authentication
- User management
- Role-based access control
"""

import hashlib
import hmac
import logging
import os
import secrets
from datetime import datetime, timedelta
from enum import Enum
from typing import Optional

from pydantic import BaseModel, EmailStr, Field

logger = logging.getLogger(__name__)

# Try to import JWT library
try:
    import jwt
    HAS_JWT = True
except ImportError:
    HAS_JWT = False
    logger.warning("PyJWT not installed. Install with: pip install PyJWT")


class UserRole(str, Enum):
    """User roles for access control."""
    ADMIN = "admin"
    USER = "user"
    READONLY = "readonly"


class User(BaseModel):
    """User model."""
    id: str
    email: EmailStr
    username: str
    role: UserRole = UserRole.USER
    is_active: bool = True
    created_at: datetime = Field(default_factory=datetime.utcnow)
    last_login: Optional[datetime] = None

    # Hashed password (never store plain text)
    password_hash: str = ""

    # API key hash for programmatic access (store hash, not plaintext)
    api_key_hash: Optional[str] = None
    api_key_prefix: Optional[str] = None  # first 8 chars for identification
    api_key_created: Optional[datetime] = None


class TokenPayload(BaseModel):
    """JWT token payload."""
    sub: str  # User ID
    email: str
    role: str
    exp: datetime
    iat: datetime = Field(default_factory=datetime.utcnow)
    type: str = "access"  # access or refresh


def _get_secret_key() -> str:
    """Get auth secret key from environment. Fail loudly if not set."""
    key = os.getenv("AUTH_SECRET_KEY")
    if not key:
        raise RuntimeError(
            "AUTH_SECRET_KEY environment variable is required. "
            "Generate one with: python -c \"import secrets; print(secrets.token_hex(32))\""
        )
    return key


class AuthConfig(BaseModel):
    """Authentication configuration."""
    secret_key: str = Field(default_factory=_get_secret_key)
    algorithm: str = "HS256"
    access_token_expire_minutes: int = 30
    refresh_token_expire_days: int = 7
    api_key_prefix: str = "alphagwas_"


class AuthManager:
    """
    Authentication manager for AlphaGWAS.

    Handles:
    - Password hashing and verification
    - JWT token creation and validation
    - API key generation and validation
    """

    def __init__(self, config: Optional[AuthConfig] = None):
        """Initialize auth manager."""
        self.config = config or AuthConfig()
        self._users: dict[str, User] = {}  # In-memory store (use database in production)
        self._api_key_index: dict[str, str] = {}  # key_hash -> user_id mapping

    # ==================== Password Management ====================

    def hash_password(self, password: str) -> str:
        """Hash password using PBKDF2."""
        salt = secrets.token_hex(16)
        hash_obj = hashlib.pbkdf2_hmac(
            'sha256',
            password.encode(),
            salt.encode(),
            100000
        )
        return f"{salt}${hash_obj.hex()}"

    def verify_password(self, password: str, password_hash: str) -> bool:
        """Verify password against hash."""
        try:
            salt, stored_hash = password_hash.split('$')
            hash_obj = hashlib.pbkdf2_hmac(
                'sha256',
                password.encode(),
                salt.encode(),
                100000
            )
            return hmac.compare_digest(hash_obj.hex(), stored_hash)
        except (ValueError, AttributeError):
            return False

    # ==================== User Management ====================

    def create_user(
        self,
        email: str,
        username: str,
        password: str,
        role: UserRole = UserRole.USER,
    ) -> User:
        """Create a new user."""
        user_id = secrets.token_hex(8)

        user = User(
            id=user_id,
            email=email,
            username=username,
            role=role,
            password_hash=self.hash_password(password),
        )

        self._users[user_id] = user
        logger.info(f"Created user: {username} ({email})")
        return user

    def get_user(self, user_id: str) -> Optional[User]:
        """Get user by ID."""
        return self._users.get(user_id)

    def get_user_by_email(self, email: str) -> Optional[User]:
        """Get user by email."""
        for user in self._users.values():
            if user.email == email:
                return user
        return None

    def get_user_by_username(self, username: str) -> Optional[User]:
        """Get user by username."""
        for user in self._users.values():
            if user.username == username:
                return user
        return None

    def authenticate_user(self, username: str, password: str) -> Optional[User]:
        """Authenticate user with username/password."""
        user = self.get_user_by_username(username)
        if not user:
            user = self.get_user_by_email(username)

        if user and user.is_active and self.verify_password(password, user.password_hash):
            user.last_login = datetime.utcnow()
            return user
        return None

    def update_password(self, user_id: str, new_password: str) -> bool:
        """Update user password."""
        user = self.get_user(user_id)
        if user:
            user.password_hash = self.hash_password(new_password)
            return True
        return False

    # ==================== JWT Tokens ====================

    def create_access_token(self, user: User) -> str:
        """Create JWT access token."""
        if not HAS_JWT:
            raise RuntimeError("PyJWT not installed")

        expire = datetime.utcnow() + timedelta(minutes=self.config.access_token_expire_minutes)

        payload = {
            "sub": user.id,
            "email": user.email,
            "role": user.role.value,
            "exp": expire,
            "iat": datetime.utcnow(),
            "type": "access",
        }

        return jwt.encode(payload, self.config.secret_key, algorithm=self.config.algorithm)

    def create_refresh_token(self, user: User) -> str:
        """Create JWT refresh token."""
        if not HAS_JWT:
            raise RuntimeError("PyJWT not installed")

        expire = datetime.utcnow() + timedelta(days=self.config.refresh_token_expire_days)

        payload = {
            "sub": user.id,
            "type": "refresh",
            "exp": expire,
            "iat": datetime.utcnow(),
        }

        return jwt.encode(payload, self.config.secret_key, algorithm=self.config.algorithm)

    def verify_token(self, token: str) -> Optional[TokenPayload]:
        """Verify and decode JWT token."""
        if not HAS_JWT:
            raise RuntimeError("PyJWT not installed")

        try:
            payload = jwt.decode(
                token,
                self.config.secret_key,
                algorithms=[self.config.algorithm]
            )
            return TokenPayload(**payload)
        except jwt.ExpiredSignatureError:
            logger.warning("Token expired")
            return None
        except jwt.InvalidTokenError as e:
            logger.warning(f"Invalid token: {e}")
            return None

    def refresh_access_token(self, refresh_token: str) -> Optional[str]:
        """Create new access token from refresh token."""
        payload = self.verify_token(refresh_token)

        if not payload or payload.type != "refresh":
            return None

        user = self.get_user(payload.sub)
        # Re-verify user is active before issuing new token
        if not user or not user.is_active:
            return None

        return self.create_access_token(user)

    # ==================== API Keys ====================

    def _hash_api_key(self, key: str) -> str:
        """Hash an API key for storage."""
        return hashlib.sha256(key.encode()).hexdigest()

    def generate_api_key(self, user: User) -> str:
        """Generate API key for user. Returns the key (only shown once)."""
        key = f"{self.config.api_key_prefix}{secrets.token_hex(24)}"
        key_hash = self._hash_api_key(key)

        # Store hash, not the key itself
        user.api_key_hash = key_hash
        user.api_key_prefix = key[:8]
        user.api_key_created = datetime.utcnow()
        self._api_key_index[key_hash] = user.id

        logger.info(f"Generated API key for user: {user.username}")
        return key

    def verify_api_key(self, api_key: str) -> Optional[User]:
        """Verify API key and return user."""
        key_hash = self._hash_api_key(api_key)
        user_id = self._api_key_index.get(key_hash)
        if user_id:
            user = self.get_user(user_id)
            if user and user.is_active and user.api_key_hash == key_hash:
                return user
        return None

    def revoke_api_key(self, user: User) -> bool:
        """Revoke user's API key."""
        if user.api_key_hash and user.api_key_hash in self._api_key_index:
            del self._api_key_index[user.api_key_hash]
            user.api_key_hash = None
            user.api_key_prefix = None
            user.api_key_created = None
            logger.info(f"Revoked API key for user: {user.username}")
            return True
        return False


# FastAPI Integration
def create_auth_routes(auth_manager: AuthManager):
    """Create FastAPI authentication routes."""
    from fastapi import APIRouter, Depends, HTTPException, status
    from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer, APIKeyHeader

    router = APIRouter(prefix="/auth", tags=["Authentication"])
    security = HTTPBearer(auto_error=False)
    api_key_header = APIKeyHeader(name="X-API-Key", auto_error=False)

    class LoginRequest(BaseModel):
        username: str
        password: str

    class TokenResponse(BaseModel):
        access_token: str
        refresh_token: str
        token_type: str = "bearer"
        expires_in: int

    class RegisterRequest(BaseModel):
        email: EmailStr
        username: str
        password: str

    async def get_current_user(
        credentials: Optional[HTTPAuthorizationCredentials] = Depends(security),
        api_key: Optional[str] = Depends(api_key_header),
    ) -> User:
        """Get current authenticated user."""
        # Try API key first
        if api_key:
            user = auth_manager.verify_api_key(api_key)
            if user:
                return user

        # Try JWT token
        if credentials:
            payload = auth_manager.verify_token(credentials.credentials)
            if payload:
                user = auth_manager.get_user(payload.sub)
                if user and user.is_active:
                    return user

        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid authentication credentials",
            headers={"WWW-Authenticate": "Bearer"},
        )

    async def require_admin(user: User = Depends(get_current_user)) -> User:
        """Require admin role."""
        if user.role != UserRole.ADMIN:
            raise HTTPException(
                status_code=status.HTTP_403_FORBIDDEN,
                detail="Admin access required",
            )
        return user

    @router.post("/register", response_model=dict)
    async def register(request: RegisterRequest):
        """Register a new user."""
        if auth_manager.get_user_by_email(request.email):
            raise HTTPException(status_code=400, detail="Email already registered")
        if auth_manager.get_user_by_username(request.username):
            raise HTTPException(status_code=400, detail="Username already taken")

        user = auth_manager.create_user(
            email=request.email,
            username=request.username,
            password=request.password,
        )

        return {"message": "User created successfully", "user_id": user.id}

    @router.post("/login", response_model=TokenResponse)
    async def login(request: LoginRequest):
        """Login and get access token."""
        user = auth_manager.authenticate_user(request.username, request.password)

        if not user:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Invalid username or password",
            )

        return TokenResponse(
            access_token=auth_manager.create_access_token(user),
            refresh_token=auth_manager.create_refresh_token(user),
            expires_in=auth_manager.config.access_token_expire_minutes * 60,
        )

    @router.post("/refresh", response_model=dict)
    async def refresh_token(refresh_token: str):
        """Refresh access token."""
        new_token = auth_manager.refresh_access_token(refresh_token)

        if not new_token:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="Invalid refresh token",
            )

        return {"access_token": new_token, "token_type": "bearer"}

    @router.get("/me", response_model=dict)
    async def get_me(user: User = Depends(get_current_user)):
        """Get current user info."""
        return {
            "id": user.id,
            "email": user.email,
            "username": user.username,
            "role": user.role.value,
            "created_at": user.created_at.isoformat(),
        }

    @router.post("/api-key", response_model=dict)
    async def create_api_key(user: User = Depends(get_current_user)):
        """Generate new API key."""
        api_key = auth_manager.generate_api_key(user)
        return {"api_key": api_key, "message": "Store this key securely - it won't be shown again"}

    @router.delete("/api-key", response_model=dict)
    async def revoke_api_key_endpoint(user: User = Depends(get_current_user)):
        """Revoke current API key."""
        if auth_manager.revoke_api_key(user):
            return {"message": "API key revoked"}
        raise HTTPException(status_code=404, detail="No API key to revoke")

    return router, get_current_user, require_admin
