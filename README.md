# Forging AES-GCM auth tags when IVs are reused

If you pass a list of ciphertexts which were produced using the same IV to `recover_auth_secret`, it will return a list of potential `(auth_key, keyblock_0)` values. These values can be passed to the `auth_tag` function to produce a valid authentication tag for an arbitrarily modified ciphertext. Remeber that this function expects a ciphertext *without* a trailing auth tag, so you probably want to strip off the last 16 bytes from a ciphertext before manipulating it and passing it to `auth_tag`.

Besides the authentication feature, GCM mode works like CTR mode or other stream ciphers, meaning flipping a bit in the ciphertext will flip the corresponding bit in the plaintext.

